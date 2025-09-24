"""Streamlit app for visualizing single-nuclide radioactive decay."""
from __future__ import annotations

import csv
import io
import math
from typing import List

import matplotlib.pyplot as plt
import streamlit as st

from src.decay import simulate_decay


def compute_analytic(N0: float, decay_constant: float, times: List[float]) -> List[float]:
    """Return analytic solution values for the provided times."""
    return [N0 * math.exp(-decay_constant * t) for t in times]


def to_csv_bytes(times: List[float], numeric: List[float], analytic: List[float]) -> bytes:
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(["시간 t (초)", "핵자수 N(t) (개)", "해석해 N₀·exp(-λt) (개)"])
    for row in zip(times, numeric, analytic):
        writer.writerow(row)
    return buffer.getvalue().encode("utf-8")


def to_preview_rows(times: List[float], numeric: List[float], analytic: List[float], limit: int = 10):
    headers = ["시간 t (초)", "핵자수 N(t) (개)", "해석해 N₀·exp(-λt) (개)"]
    rows = []
    for t, n, a in list(zip(times, numeric, analytic))[:limit]:
        rows.append(dict(zip(headers, (t, n, a))))
    return rows


st.set_page_config(page_title="단일 핵종 붕괴 시각화", layout="wide")
st.title("단일 핵종 방사성 붕괴 시각화")

with st.sidebar:
    st.header("시뮬레이션 설정")
    T12 = st.slider("반감기 T₁/₂ (초)", min_value=0.1, max_value=100.0, value=10.0, step=0.1)
    N0 = st.slider("초기 핵자수 N₀ (개)", min_value=10, max_value=10000, value=1000, step=10)
    T_total = st.slider("총 시간 T (초)", min_value=1.0, max_value=300.0, value=60.0, step=1.0)
    dt = st.slider("시간 간격 Δt (초)", min_value=0.1, max_value=10.0, value=0.5, step=0.1)

lambda_decay = math.log(2.0) / T12

times, numeric_N = simulate_decay(T12=T12, N0=N0, dt=dt, T_total=T_total)
analytic_N = compute_analytic(N0=N0, decay_constant=lambda_decay, times=times)

csv_bytes = to_csv_bytes(times, numeric_N, analytic_N)
preview_rows = to_preview_rows(times, numeric_N, analytic_N)

with st.sidebar:
    st.download_button(
        label="CSV로 저장",
        data=csv_bytes,
        file_name="decay_simulation.csv",
        mime="text/csv",
    )

col_plot, col_table = st.columns([2, 1])

with col_plot:
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times, numeric_N, label="수치해 N(t)", color="tab:blue")
    ax.plot(times, analytic_N, label="해석해 N₀·exp(-λt)", color="tab:orange", linestyle="--")
    ax.set_xlabel("시간 t (초)")
    ax.set_ylabel("핵자수 N (개)")
    ax.set_title("방사성 붕괴 곡선")
    ax.legend()
    ax.grid(True, linestyle=":", linewidth=0.5)
    st.pyplot(fig)

with col_table:
    st.subheader("데이터 미리보기")
    st.table(preview_rows)

st.caption(
    "λ = ln(2) / T₁/₂, N(t) = N₀ · exp(-λt). 시간 단위: 초, 핵자수 단위: 개."
)
