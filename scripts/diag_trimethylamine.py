#!/usr/bin/env python
"""扫描 Py-GC-MS 处理样品的 TMAH 副产物（三甲胺，基峰 m/z 58）峰。

严谨性背景：TMAH（四甲基氢氧化铵）热化学分解产物会被 NIST 误命名为各种
含 N 化合物并落入 Other_N -> Microbial，虚假抬高微生物信号与 R_MP。
基于 *名称* 的剔除（pipeline TMAH_ARTIFACTS）可能漏掉误命名的试剂峰；
本脚本基于 *谱* 判定：基峰 m/z 58±0.5 且 59/42 碎片成比例 -> 三甲胺特征，
无论 NIST 给它什么名字都标记为 TMAH 试剂峰。

用法:
    python diag_trimethylamine.py --qgd <QGD_dir> --txt <TXT_dir> \
        --sample_map sample_map.json [--rt_tolerance 0.35] [--tic_min 50000]

    <QGD_dir> 内含 {sid}.qgd；<TXT_dir> 内含 {sid}.txt（NIST 导出）。
    sample_map.json: {'5':'CK','6':'BC7.5',...}

输出: 控制台逐峰报告（QGD_RT / TIC / frac58 / 谱 / TXT 匹配名 / SI / 分类），
      以及可选的 --output JSON 清单（推荐存为 corrections.json 输入给 pipeline）。
"""
import os, sys, json, argparse

# 使脚本可从 skill 目录任意位置运行
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qgd_reader import QGDFile
from pipeline import parse_txt, classify_compound, load_shahriar_library


def tma_score(peaks):
    """返回 (is_tma, frac58) — 基峰 m/z 58 且 59/42 碎片成比例。"""
    # 排除常见背景峰（N2 m/z28、CO2 m/z44/45）
    spec = [p for p in peaks if round(p[0], 1) not in {28.0, 43.9, 44.0, 44.1, 45.0}]
    if not spec:
        return False, 0.0
    top = sorted(spec, key=lambda x: -x[1])
    base_mz, base_i = top[0]
    if not (57.5 <= base_mz <= 58.5):
        return False, 0.0
    total = sum(i for _, i in spec) or 1
    frac58 = base_i / total
    i59 = next((i for mz, i in top if 58.5 <= mz <= 59.5), 0)
    i42 = next((i for mz, i in top if 41.5 <= mz <= 42.5), 0)
    ok59 = i59 / base_i > 0.10
    ok42 = i42 / base_i > 0.15
    return (frac58 > 0.12 and ok59 and ok42), frac58


def scan(qgd_path, txt_path, label, rt_tol, tic_min, lib):
    """扫描单个样品：QGD TIC 强峰 -> 峰顶谱 -> TMA 判定 -> TXT 匹配。"""
    q = QGDFile(qgd_path)
    tic = q.get_tic_data()  # [(rt, tic)]
    txt_peaks = parse_txt(txt_path) if os.path.exists(txt_path) else []

    strong = sorted(tic, key=lambda x: -x[1])[:60]
    hits = []
    seen = set()
    for rt, tic_i in strong:
        local = [(r, t) for r, t in tic if abs(r - rt) < 0.15]
        if not local:
            continue
        mx = max(t for _, t in local)
        if tic_i < mx * 0.9:  # 非局部极值
            continue
        key = round(rt, 2)
        if key in seen:
            continue
        seen.add(key)
        if tic_i < tic_min:
            continue
        s = q.get_spectrum_at_rt(rt, tolerance=0.06)
        if not s:
            continue
        is_tma, frac58 = tma_score(s["peaks"])
        if is_tma:
            hits.append((rt, tic_i, frac58, s["peaks"]))

    print(f"\n--- {label} ---")
    if not hits:
        print("  无 TMA（基峰58）特征峰")
        q.close()
        return []

    results = []
    for rt_q, tic_i, frac, spec in hits:
        # TXT 匹配：时间基偏移窗口
        best, best_d = None, rt_tol
        for p in txt_peaks:
            d = abs(p["rt"] - rt_q)
            if d < best_d and p["conc"] > 0.3:
                best_d = d
                best = p
        mz_top3 = ", ".join(f"{m:.0f}({i/1e6:.1f}M)" for m, i in sorted(spec, key=lambda x: -x[1])[:3])
        entry = {"qgd_rt": round(rt_q, 2), "tic": tic_i, "frac58": round(frac, 2),
                 "spectrum": mz_top3}
        if best is None:
            entry["txt_rt"] = None
            entry["name"] = None
            entry["status"] = "QGD 无 TXT 显著匹配"
            print(f"  QGD_RT={rt_q:.2f} TIC={tic_i:,.0f} (frac58={frac:.2f}) — TXT 无匹配显著峰")
        else:
            cls = classify_compound(best["name"], lib)
            entry.update({"txt_rt": round(best["rt"], 3), "conc": round(best["conc"], 2),
                          "si": best["si"], "class": cls, "name": best["name"],
                          "status": "TMAH_REAGENT_PEAK"})
            print(f"  QGD_RT={rt_q:.2f} TIC={tic_i:,.0f} frac58={frac:.2f} 谱={mz_top3}")
            print(f"     -> TXT_RT={best['rt']:.3f} Conc={best['conc']:.2f} SI={best['si']} Class={cls}")
            print(f"       Name={best['name'][:80]}")
        results.append(entry)
    q.close()
    return results


def main():
    ap = argparse.ArgumentParser(description="TMAH 试剂峰谱检（基峰 m/z 58）")
    ap.add_argument("--qgd", required=True, help="QGD 文件目录")
    ap.add_argument("--txt", required=True, help="NIST TXT 目录")
    ap.add_argument("--sample_map", help="JSON: {'5':'CK',...}")
    ap.add_argument("--features", help="可选：Stage1 features_clean.csv（含 peak_id，用于关联决策）")
    ap.add_argument("--rt_tolerance", type=float, default=0.35,
                    help="TXT/QGD 时间基偏移容差 (min)")
    ap.add_argument("--tic_min", type=float, default=50000.0,
                    help="TIC 峰高下限")
    ap.add_argument("--output", help="输出目录（写入 tmah_decisions.csv 与 JSON 清单）")
    args = ap.parse_args()

    lib = load_shahriar_library()
    sample_map = {"5": "CK", "6": "BC7.5", "7": "BC15", "8": "BC30"}
    if args.sample_map and os.path.exists(args.sample_map):
        with open(args.sample_map) as f:
            sample_map = json.load(f)
            sample_map.pop("_notes", None)

    # 预加载 features_clean.csv（若提供）用于 peak_id 关联
    peak_by_key = {}  # (sample_id, rt) -> peak_id
    if args.features and os.path.exists(args.features):
        import csv as _csv
        with open(args.features, encoding="utf-8-sig") as f:
            for row in _csv.DictReader(f):
                try:
                    rt = float(row.get("rt_min", ""))
                except (TypeError, ValueError):
                    continue
                peak_by_key[(row.get("sample_id", ""), round(rt, 3))] = row.get("peak_id", "")
        print(f"Loaded {len(peak_by_key)} peaks from features_clean.csv")

    all_hits = {}
    decision_rows = []
    for sid in sorted(sample_map.keys()):
        trt = sample_map[sid]
        qgd = os.path.join(args.qgd, f"{sid}.qgd")
        if not os.path.exists(qgd):
            qgd = os.path.join(args.qgd, f"{sid}.QGD")
        if not os.path.exists(qgd):
            print(f"\n  {trt}: QGD 不存在（跳过）")
            continue
        txt = os.path.join(args.txt, f"{sid}.txt")
        if not os.path.exists(txt):
            txt = os.path.join(args.txt, f"{sid}.TXT")
        hits = scan(qgd, txt, trt, args.rt_tolerance, args.tic_min, lib)
        for h in hits:
            if h.get("txt_rt") is None:
                continue
            rt_key = round(h["txt_rt"], 3)
            peak_id = peak_by_key.get((trt, rt_key), "")
            decision = "EXCLUDE"
            reason = "TMAH_REAGENT_PATTERN"
            # 仅当谱以 m/z58/59/42 为主时 EXCLUDE；否则 REVIEW
            if h.get("frac58", 0) < 0.20:
                decision, reason = "REVIEW", "WEAK_MZ58_SIGNAL"
            all_hits.setdefault(trt, {})[str(rt_key)] = h["name"]
            decision_rows.append({
                "peak_id": peak_id,
                "sample_id": trt,
                "rt_min": rt_key,
                "compound_name": h.get("name", ""),
                "mz58_fraction": h.get("frac58", ""),
                "tic": h.get("tic", ""),
                "decision": decision,
                "reason": reason,
            })

    if args.output:
        os.makedirs(args.output, exist_ok=True)
        # canonical Stage 4 output: tmah_decisions.csv
        import csv as _csv
        dec_path = os.path.join(args.output, "tmah_decisions.csv")
        if decision_rows:
            with open(dec_path, "w", newline="", encoding="utf-8-sig") as f:
                w = _csv.DictWriter(f, fieldnames=list(decision_rows[0].keys()))
                w.writeheader()
                for r in decision_rows:
                    w.writerow(r)
            print(f"\nTMAH 决策表已写入: {dec_path}")
        else:
            # 无检出时也写空表（带表头），保证数据契约存在
            with open(dec_path, "w", newline="", encoding="utf-8-sig") as f:
                w = _csv.DictWriter(f, fieldnames=[
                    "peak_id", "sample_id", "rt_min", "compound_name",
                    "mz58_fraction", "tic", "decision", "reason"])
                w.writeheader()
            print(f"\nTMAH 决策表（空）已写入: {dec_path}")
        json_path = os.path.join(args.output, "tmah_hits.json")
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(all_hits, f, ensure_ascii=False, indent=2)
        print(f"TMAH 谱检清单已写入: {json_path}")
        print("注意: 请勿将剔除语义写入 corrections.json（该文件仅用于改名）。"
              "Stage 5 apply_final.py 会读取 tmah_decisions.csv 执行剔除。")


if __name__ == "__main__":
    main()
