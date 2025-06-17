import os
import re
import numpy as np
import meshio
import matplotlib.pyplot as plt

# ─ 사용자 설정 ─
vtk_dir = "/Users/eunchan/workspace/HW4/outputs"
varname = "velocity"  # C 코드에서 print_vtk_vector() 호출 시 사용한 varname
Re = 500  # 동일하게 _CONST_Re 값
step_interval = 10  # 몇 스텝 간격으로 뽑을지

# 파일명 패턴: e.g. velocity_Re100_t0000.vtk
pattern = re.compile(rf"{varname}_Re{Re}_t(\d{{4}})\.vtk")

iterations = []
residuals_u = []
residuals_v = []

prev_u = None
prev_v = None

for fname in sorted(os.listdir(vtk_dir)):
    m = pattern.match(fname)
    if not m:
        continue
    step = int(m.group(1))
    if step % step_interval != 0:
        continue

    path = os.path.join(vtk_dir, fname)
    mesh = meshio.read(path)
    # VTK VECTORS <varname> double → shape (N_points, 3)
    vec = mesh.point_data[varname]
    u = vec[:, 0]
    v = vec[:, 1]

    if prev_u is not None:
        # 이전 스텝과 비교한 차이 절댓값의 합 (L1-norm)
        res_u = np.sum(np.abs(u - prev_u)) / np.sum(np.abs(u))
        res_v = np.sum(np.abs(v - prev_v)) / np.sum(np.abs(v))

        iterations.append(step)
        residuals_u.append(res_u)
        residuals_v.append(res_v)

    prev_u = u.copy()
    prev_v = v.copy()

# ─ 결과 플롯 ─
plt.figure(figsize=(8, 5))
plt.plot(iterations, residuals_u, "-", label="u")
plt.plot(iterations, residuals_v, "--", label="v")
plt.xlabel("Iteration")
plt.ylabel("Residuals (L1-norm)")
plt.title(f"{varname} Residuals (Re={Re}, every {step_interval} steps)")
plt.yscale("log")
plt.legend()
plt.grid(True)
plt.tight_layout()
# Save the figure to file
plt.savefig(f"{varname}_res_Re{Re}_step{step_interval}.png", dpi=300)
# plt.show()
