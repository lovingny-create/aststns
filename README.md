# Astrophysics/Climate Simulation Toolkit

이 저장소에는 식쌍성 광도곡선, 지구 공전·일사량(밀란코비치 주기 포함), 천정 성도 퀴즈 등 여러 교육용 시뮬레이터와 스크립트가 포함되어 있습니다.

## 구성
- `index.html` / `style.css` / `script.js`: Plotly.js 그래프와 슬라이더 UI를 갖춘 실시간 식쌍성(eclipsing binary) 광도곡선 웹 시뮬레이터
- `binary_lightcurve.py`: 별, 식쌍성 궤도, 광도 곡선 계산 로직을 담은 모듈
- `simulate.py`: 명령행에서 식쌍성 시뮬레이션을 실행하고 CSV·PNG 그래프로 광도 곡선을 저장하는 스크립트
- `app.py`: Streamlit 기반으로 슬라이더를 제공하는 실시간 식쌍성 데스크톱 웹앱
- `solar_insolation_app.py`: 지구 공전 궤도, 세차각(ω), 이심률(e), 축경사(ε)가 바뀔 때의 적위·남중고도·위도별 일사량을 확인하는 Streamlit 시뮬레이터 (밀란코비치 주기 탐구용)
- `missing_star_app.py`: 임의의 별을 숨긴 천정(zenith) 성도 문제·정답 이미지를 생성하는 Streamlit 앱

## 빠른 시작

### 식쌍성 광도곡선 (웹)
1. 저장소 루트에서 `index.html`을 브라우저로 직접 열거나 정적 서버로 제공합니다.
2. 반지름(r1, r2), 표면온도(T1, T2), 질량(M1, M2), 이심률(e), 경사각(i), 궤도주기(P) 슬라이더를 조정해 Plotly.js 그래프를 실시간으로 확인합니다.

### 식쌍성 광도곡선 (Streamlit)
1. 의존성 설치: `pip install streamlit plotly`
2. 앱 실행: `streamlit run app.py`
3. 사이드바에서 반지름, 온도, 경사각, 이심률, 질량, 장반경, 샘플 개수를 조정하면 그래프와 주기가 즉시 갱신됩니다.

### 지구 공전·밀란코비치 주기 시뮬레이터 (Streamlit)
1. 의존성 설치: `pip install -r requirements.txt`
2. 실행: `streamlit run solar_insolation_app.py`
3. 달·일 또는 N일차를 선택하고 이심률, 세차각(ω), 축경사(ε), 위도를 조정해 궤도 그림, 적위, 위도별 일사량 곡선, 남중고도를 확인합니다.

### 별 숨기기 성도 생성기 (Streamlit)
1. 의존성 설치: `pip install -r requirements.txt` (추가로 `starplot` 필요)
2. 실행: `streamlit run missing_star_app.py`
3. 날짜/시간, 관측 위치(위도·경도), 밝기 등급 한계를 입력 후 “👉 성도 생성하기”를 클릭하면 문제·정답 성도 PNG와 숨겨진 별 목록이 표시됩니다.

### CLI로 식쌍성 CSV 저장
```bash
python simulate.py --mass1 1.989e30 --mass2 1.5e30 \
  --radius1 6.957e8 --radius2 5e8 \
  --temp1 6000 --temp2 5000 \
  --semi-major-axis 2.5e10 --eccentricity 0.1 \
  --inclination $(python - <<'PY'
import math
print(math.radians(87))
PY
) \
  --phases 400 --output light_curve.csv
```

실행 후 `light_curve.csv`에는 위상(0~1)과 정규화된 상대 광도가 저장되고, 동일한 이름의 PNG(`light_curve.png`)로 광도 곡선 플롯이 생성됩니다. `--figure`로 PNG 경로를 바꿀 수 있으며, 표준 출력에는 궤도 주기(일 단위)도 표시됩니다. 시스템에 matplotlib이 없다면 같은 이름의 SVG로 대체 저장합니다.

## 모델 가정 및 한계
- 별은 완전한 구 형태이며 림브 다크닝, 타이달 변형, 반사 효과 등을 고려하지 않습니다.
- 궤도는 점질량 두 개의 케플러 운동으로 모델링하며, 분광선 속도나 시간 지연 효과는 포함하지 않습니다.
- 광도는 각 별의 표면 밝기를 스테판-볼츠만 법칙으로 계산한 후 식 현상 시 기하학적 가려짐에 따른 면적 감소만을 반영합니다.
- 일사량 시뮬레이터는 일 평균값을 사용하며 대기 효과나 구름, 지표 알베도는 고려하지 않습니다.

## 개발/테스트
```bash
python -m py_compile binary_lightcurve.py simulate.py solar_insolation_app.py
```

## 고유운동·연주시차 시뮬레이터 미리보기
- 저장소 루트에서 `python -m http.server 8000`를 실행한 뒤 브라우저에서 `http://localhost:8000/proper_motion_parallax.html`로 접속하면 전체 UI를 바로 확인할 수 있습니다.
- 파일을 직접 더블클릭해 열어도 동일하게 작동하지만, 로컬 서버를 거치면 외부 리소스 접근 문제 없이 일정하게 표시됩니다.
- 기본 제공되는 슬라이더 값(연주시차 100 mas, β=30° 등)으로 로드되며, 좌측 카드에서 파라미터를 조정하면 오른쪽 궤적 플롯이 즉시 갱신됩니다. 애니메이션 하단 슬라이더로 특정 시점을 재생하거나 CSV 스케줄 파일을 불러와 관측 시각을 테스트할 수 있습니다.
