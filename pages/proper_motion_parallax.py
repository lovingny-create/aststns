"""Streamlit wrapper for the proper motion + annual parallax simulator."""
from __future__ import annotations

from pathlib import Path

import streamlit as st


def load_html() -> str:
    """Load the standalone HTML file that renders the simulator.

    Keeping the HTML as a separate asset makes it easier to iterate on the
    visualization without changing the Streamlit wrapper.
    """

    html_path = Path(__file__).resolve().parent.parent / "proper_motion_parallax.html"
    return html_path.read_text(encoding="utf-8")


def main() -> None:
    st.set_page_config(page_title="고유운동 + 연주시차 시뮬레이터", layout="wide")
    st.title("고유운동 + 연주시차 시뮬레이터")

    st.markdown(
        """
        기존 단일 HTML로 만든 고유운동·연주시차 시뮬레이터를 Streamlit 페이지로 감싸
        앱으로 배포할 수 있게 했습니다. 아래 포함된 뷰는 원본 HTML/JS를 그대로 사용하므로,
        슬라이더와 라디오 버튼을 조정하며 연주시차 타원, 고유운동, 샘플 스케줄을 즉시 확인할 수 있습니다.
        """
    )

    st.components.v1.html(load_html(), height=1100, scrolling=True)


if __name__ == "__main__":
    main()
