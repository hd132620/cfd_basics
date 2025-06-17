# fc_ab2_cn_adi

2D Navier-Stokes 및 Poisson 방정식 해석용 CFD (Lid-driven cavity) 연습코드
![alt text](image.png)
---

## 수치 해법 요약

- 공간 이산화: 2nd order central difference
- 점성항 시간 적분: Crank-Nicolson (2nd, implicit)
- 대류항 시간 적분: Adams-Bashforth (2nd, explicit)
- 압력-속도 결합: Fractional Step Method (Face-centered pressure)
- 선형 시스템 풀이: ADI(Alternating Direction Implicit)-TDMA

---

## 주요 파일

- `main.c` : 실행 및 전체 흐름
- `navier_stokes.c/h` : Navier-Stokes 해법
- `poisson.c/h` : Poisson_ADI / CN_ADI 해법
- `setting.c/h` : 격자/초기/경계 조건
- `outputs/` : 결과 저장
- `postprocessing/plot_res.py` : 결과 시각화

---

## 실행

```sh
make
./main
cd postprocessing
python3 plot_res.py
```
