"""
Shimadzu QGD File Reader for Py-GC-MS
=====================================
Extracts TIC, mass spectra, and metadata from Shimadzu GCMSsolution .qgd files.

QGD format (OLE2 Compound Document):
  - GCMS Raw Data/Retention Time: int32[scan_count] in milliseconds
  - GCMS Raw Data/Spectrum Index: int32[scan_count] byte offsets into MS Raw Data
  - GCMS Raw Data/TIC Data: int64[scan_count] TIC intensities
  - GCMS Raw Data/MS Raw Data: concatenated mass spectra

Each mass spectrum block (32-byte header):
  [0:4]   scan_number   (int32)
  [4:8]   rt_ms          (int32, /60000 for minutes)
  [8:20]  flags/skip
  [20:22] n_bytes        (int16) - bytes per intensity value (1-5)
  [22:24] nval           (int16) - number of m/z peaks
  [24:32] skip
  Then: nval x (2-byte m/z_raw [int16, /20 for actual m/z] + n_bytes intensity)

Based on chromConverter by Ethan Bass (CC-BY).
"""
import olefile
import struct
import os
import csv


class QGDFile:
    """Parser for Shimadzu GCMSsolution .qgd files."""

    def __init__(self, filepath):
        self.filepath = filepath
        self.ole = olefile.OleFileIO(filepath)
        self._read_scan_metadata()
        self._read_file_properties()

    def _read_scan_metadata(self):
        """Parse retention times, spectrum indices, and TIC."""
        rt_bytes = self.ole.openstream("GCMS Raw Data/Retention Time").read()
        self.scan_count = len(rt_bytes) // 4
        self.rt_min = [
            struct.unpack_from("<i", rt_bytes, i * 4)[0] / 60000.0
            for i in range(self.scan_count)
        ]
        idx_bytes = self.ole.openstream("GCMS Raw Data/Spectrum Index").read()
        self.spec_offsets = [
            struct.unpack_from("<i", idx_bytes, i * 4)[0]
            for i in range(self.scan_count)
        ]
        tic_bytes = self.ole.openstream("GCMS Raw Data/TIC Data").read()
        self.tic = [
            struct.unpack_from("<q", tic_bytes, i * 8)[0]
            for i in range(self.scan_count)
        ]
        self.ms_raw = self.ole.openstream("GCMS Raw Data/MS Raw Data").read()
        self._spec_end = self.spec_offsets + [len(self.ms_raw)]

    def _read_file_properties(self):
        """Extract metadata from File Property stream."""
        self.metadata = {}
        try:
            fp = self.ole.openstream("File Property").read()
            for key, offset in [("sample_name", 204), ("operator", 300), ("datafile", 580)]:
                try:
                    end = fp.index(b"\x00", offset)
                    self.metadata[key] = fp[offset:end].decode("utf-16-le", errors="replace")
                except:
                    pass
        except:
            pass

    def get_spectrum(self, scan_index):
        """Parse a single mass spectrum by scan index.

        Returns dict with: scan, rt, n_peaks, peaks=[(mz, intensity), ...]
        """
        if scan_index >= self.scan_count:
            return None
        start = self._spec_end[scan_index]
        end = self._spec_end[scan_index + 1]
        data = self.ms_raw[start:end]
        if len(data) < 32:
            return None

        rt_ms = struct.unpack_from("<i", data, 4)[0]
        n_bytes = struct.unpack_from("<h", data, 20)[0]
        nval = struct.unpack_from("<h", data, 22)[0]

        # Validate and recover
        expected = 32 + nval * (2 + n_bytes)
        if expected != len(data):
            for try_nb in range(1, 6):
                if (len(data) - 32) % (2 + try_nb) == 0:
                    n_bytes = try_nb
                    nval = (len(data) - 32) // (2 + try_nb)
                    break

        peaks = []
        offset = 32
        add_byte = n_bytes % 2 == 1
        nb_read = n_bytes + 1 if add_byte else n_bytes

        for _ in range(nval):
            if offset + 2 + n_bytes > len(data):
                break
            mz = struct.unpack_from("<h", data, offset)[0] / 20.0
            offset += 2
            raw = data[offset : offset + n_bytes]
            if add_byte:
                raw = raw + b"\x00"
            fmt = {1: "<B", 2: "<H", 4: "<i", 6: "<q", 8: "<q"}.get(nb_read, "<q")
            raw_pad = raw[:nb_read] if len(raw) >= nb_read else raw + b"\x00" * (nb_read - len(raw))
            intensity = struct.unpack(fmt, raw_pad)[0]
            if n_bytes == 4:
                intensity &= 0x7FFFFFFF
            offset += n_bytes
            if mz > 0:
                peaks.append((round(mz, 1), intensity))

        return {"scan": scan_index, "rt": rt_ms / 60000.0, "n_peaks": len(peaks), "peaks": peaks}

    def get_spectrum_at_rt(self, target_rt, tolerance=0.03):
        """Get spectrum closest to a target retention time (in minutes)."""
        best = min(range(self.scan_count), key=lambda i: abs(self.rt_min[i] - target_rt))
        if abs(self.rt_min[best] - target_rt) > tolerance:
            return None
        return self.get_spectrum(best)

    def get_tic_data(self):
        """Return TIC as list of (rt_min, intensity) tuples."""
        return list(zip(self.rt_min, self.tic))

    def export_tic_csv(self, path):
        """Export TIC to CSV."""
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["RT_min", "TIC"])
            for rt, tic in zip(self.rt_min, self.tic):
                w.writerow([rt, tic])

    def close(self):
        self.ole.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else input("QGD file path: ")
    with QGDFile(path) as qgd:
        print(f"Scans: {qgd.scan_count}, RT: {qgd.rt_min[0]:.3f}-{qgd.rt_min[-1]:.3f} min")
        for k, v in qgd.metadata.items():
            print(f"  {k}: {v}")
        if len(sys.argv) > 2:
            spec = qgd.get_spectrum_at_rt(float(sys.argv[2]))
            if spec:
                top5 = sorted(spec["peaks"], key=lambda x: x[1], reverse=True)[:5]
                for mz, i in top5:
                    print(f"  m/z {mz:.1f}: {i:,}")
