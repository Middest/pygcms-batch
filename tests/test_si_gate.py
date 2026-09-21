#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SI 硬门槛回归测试（pipeline.py）

背景：`filter_peaks` 原本写作
    if p["si"] > 0 and p["si"] < si_threshold:
两个缺陷：
  1. `p["si"] > 0` 使 SI=0/缺失（库检索无命中）的峰**永不因 SI 被剔除** ——
     未鉴定峰被当作合格峰放进分析；
  2. 只支持 '>=' 语义，无法表达严格 '>'；而 '>=' 与 '>' 在本项目数据集上
     相差 53 个峰（3.34%），会静默改变论文里的每一个百分比。

运行: python tests/test_si_gate.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("pipeline", ROOT / "scripts" / "pipeline.py")
pl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pl)

_n = 0


def ok(cond, msg):
    global _n
    _n += 1
    if not cond:
        raise AssertionError(f"FAIL: {msg}")
    print(f"  ok {_n:02d} {msg}")


def pk(si, name="Compound", area=100.0):
    return {"si": si, "name": name, "area": area, "conc": 1.0, "rt": 5.0}


def filt(peaks, thr, op):
    return pl.filter_peaks(peaks, thr, op, remove_contaminants=False, remove_tmah=False)


def test_si_passes():
    print("si_passes:")
    ok(pl.DEFAULT_SI_THRESHOLD == 80, "默认门槛 = 80（本项目口径）")
    ok(set(pl.SI_OPERATORS) == {">=", ">"}, "仅允许 '>=' 与 '>'")
    ok(pl.si_passes(80, 80, ">="), "'>=' 保留 SI=80")
    ok(not pl.si_passes(80, 80, ">"), "'>' 剔除 SI=80（实测差 53 峰 / 3.34%）")
    ok(pl.si_passes(81, 80, ">"), "'>' 保留 SI=81")
    ok(not pl.si_passes(79, 80, ">="), "79 不达标")
    ok(pl.si_passes(50, None), "threshold=None → 不设门槛，全部通过")
    ok(not pl.si_passes(0, 80, ">="), "SI=0（库检索无命中）不达标 —— 修复点 1")
    ok(not pl.si_passes(None, 80, ">="), "SI=None 不达标 —— 修复点 1")
    ok(not pl.si_passes("", 80, ">="), "空 SI 不达标")
    ok(not pl.si_passes("abc", 80, ">="), "非数值 SI 不达标")
    try:
        pl.si_passes(90, 80, "==")
        ok(False, "非法算子应抛错")
    except ValueError:
        ok(True, "非法算子抛 ValueError")


def test_filter_peaks_gate():
    print("filter_peaks SI 硬门槛:")
    peaks = [pk(95), pk(90), pk(80), pk(79), pk(70), pk(0, ""), pk(None, "")]

    kept, rem = filt(peaks, 80, ">=")
    ok(len(kept) == 3, "SI>=80 保留 3 个（95/90/80）")
    ok(len(rem) == 4, "剔除 4 个（79/70/SI=0/SI=None）")
    reasons = " | ".join(r["filter_reason"] for r in rem)
    ok("SI=0/NA" in reasons, f"SI=0/None 被明确标注为未鉴定: {reasons}")

    kept, rem = filt(peaks, 80, ">")
    ok(len(kept) == 2, "严格 SI>80 保留 2 个（95/90）")

    kept, rem = filt(peaks, None, ">=")
    ok(len(kept) == 7, "不设门槛时全部保留（含 SI=0/None）")

    kept, rem = filt(peaks, 70, ">=")
    ok(len(kept) == 5, "SI>=70 保留 5 个")

    # 与旧行为对比：旧式 `si>0 and si<thr` 会放过 SI=0/None
    old_kept = [p for p in peaks if not (p["si"] and 0 < p["si"] < 80)]
    ok(len(old_kept) == 5, "旧逻辑会放过 5 个（含 SI=0/None）—— 本测试锁定已修复")
    ok(len(kept) < len(old_kept) + 2,
       "新逻辑不再让未鉴定峰通过")


def test_reason_labels():
    print("剔除理由与口径标签:")
    ok(pl.si_gate_label(80, ">=") == "SI>=80", "标签 SI>=80")
    ok(pl.si_gate_label(80, ">") == "SI>80", "标签 SI>80")
    ok(pl.si_gate_label(None) == "no SI gate", "标签 no SI gate")
    _, rem = filt([pk(50)], 80, ">=")
    ok("SI>=80" in rem[0]["filter_reason"], f"理由含实际口径: {rem[0]['filter_reason']}")
    _, rem = filt([pk(50)], 80, ">")
    ok("SI>80" in rem[0]["filter_reason"], f"算子在理由中可区分: {rem[0]['filter_reason']}")


def test_other_filters_unaffected():
    print("其余过滤逻辑未被破坏:")
    p = pk(95, "Bis(2-ethylhexyl) phthalate")
    kept, rem = pl.filter_peaks([p], 80, ">=")
    ok(len(rem) == 1, "污染物仍被剔除")
    # 用名单里真实存在的条目（pipeline.TMAH_ARTIFACTS）
    p2 = pk(95, "Methylamine, N,N-dimethyl-")
    kept, rem = pl.filter_peaks([p2], 80, ">=", remove_tmah=True)
    ok(len(rem) == 1, f"TMAH 副产物仍被剔除: {rem[0]['filter_reason'] if rem else '未剔除'}")
    ok(pl.TMAH_ARTIFACTS and all(isinstance(k, str) for k in pl.TMAH_ARTIFACTS),
       "TMAH_ARTIFACTS 名单存在")
    kept, rem = pl.filter_peaks([p2], 80, ">=", remove_tmah=False)
    ok(len(kept) == 1, "remove_tmah=False 时保留")
    p3 = pk(95)
    p3["area"] = 0.0
    kept, rem = pl.filter_peaks([p3], 80, ">=", remove_contaminants=False, remove_tmah=False)
    ok(len(rem) == 1 and "zero_area" in rem[0]["filter_reason"], "零面积仍被剔除")
    # SI 不达标 + 副产物 同时命中时，两个理由都要记录
    p4 = pk(50, "Methenamine")
    kept, rem = pl.filter_peaks([p4], 80, ">=", remove_tmah=True)
    ok(len(rem) == 1, "低 SI 的副产物被剔除")
    ok("SI=50" in rem[0]["filter_reason"] and "artifact" in rem[0]["filter_reason"],
       f"同时记录 SI 与副产物两个理由: {rem[0]['filter_reason']}")


def test_contract():
    print("接口契约（工作流 G1 依赖）:")
    import inspect
    for fn in (pl.run_pipeline, pl.filter_peaks):
        params = inspect.signature(fn).parameters
        ok("si_threshold" in params, f"{fn.__name__}() 暴露 si_threshold")
    ok("si_operator" in inspect.signature(pl.filter_peaks).parameters,
       "filter_peaks() 暴露 si_operator")
    ok("si_operator" in inspect.signature(pl.run_pipeline).parameters,
       "run_pipeline() 暴露 si_operator")


if __name__ == "__main__":
    print("=" * 60)
    print("pipeline.py SI 硬门槛测试")
    print("=" * 60)
    test_si_passes()
    test_filter_peaks_gate()
    test_reason_labels()
    test_other_filters_unaffected()
    test_contract()
    print("=" * 60)
    print(f"全部通过: {_n} 项断言")
