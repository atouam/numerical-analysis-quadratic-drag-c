# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# %%
cases = [
    (25 / 50, "50"),
    (25 / 100, "100"),
    (25 / 400, "400"),
]


# %%
analytic = np.loadtxt("data/vel_analytic.dat", delimiter=";")
Xa = analytic[:, 0]
Ya = analytic[:, 1]

fig = plt.figure(figsize=(24, 6))
fig.suptitle(
    "Time evolution of velocity: comparison between analytical solution, explicit Euler, and RK2 schemes for differents time steps",
    fontsize=16,
)
gs = GridSpec(1, 3, figure=fig)

for i, (dt, tag) in enumerate(cases):
    # Load Euler
    euler = np.loadtxt(f"data/vel_euler_{tag}.dat", delimiter=";")
    Xe, Ye = euler[:, 0], euler[:, 1]

    # Load RK2
    rk2 = np.loadtxt(f"data/vel_rk2_{tag}.dat", delimiter=";")
    Xr, Yr = rk2[:, 0], rk2[:, 1]

    ax = fig.add_subplot(gs[0, i])
    ax.plot(Xe, Ye, lw=2, marker=".", label="Euler")
    ax.plot(Xr, Yr, lw=2, marker=".", label="RK2")
    ax.plot(Xa, Ya, lw=2, label="Analytic")

    ax.set_title(f"dt = {dt:.4f} s")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("v (m/s)")
    ax.grid()
    ax.legend()

plt.tight_layout()
plt.savefig("figures/Velocity_schemes_comparison.png")
plt.show()


# %%
analytic = np.loadtxt("data/pos_analytic.dat", delimiter=";")
Xa = analytic[:, 0]
Ya = analytic[:, 1]

fig = plt.figure(figsize=(24, 6))
fig.suptitle(
    "Time evolution of position: comparison between analytical solution, explicit Euler, and RK2 schemes for differents time steps",
    fontsize=16,
)
gs = GridSpec(1, 3, figure=fig)

for i, (dt, tag) in enumerate(cases):
    # Load Euler
    euler = np.loadtxt(f"data/pos_euler_{tag}.dat", delimiter=";")
    Xe, Ye = euler[:, 0], euler[:, 1]

    # Load RK2
    rk2 = np.loadtxt(f"data/pos_rk2_{tag}.dat", delimiter=";")
    Xr, Yr = rk2[:, 0], rk2[:, 1]

    ax = fig.add_subplot(gs[0, i])
    ax.plot(Xe, Ye, lw=2, marker=".", label="Euler")
    ax.plot(Xr, Yr, lw=2, marker=".", label="RK2")
    ax.plot(Xa, Ya, lw=2, label="Analytic")

    ax.set_title(f"dt = {dt:.4f} s")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("v (m/s)")
    ax.grid()
    ax.legend()

plt.tight_layout()
plt.savefig("figures/Position_schemes_comparison.png")
plt.show()


# %%
analytic = np.loadtxt("data/accel_analytic.dat", delimiter=";")
Xa = analytic[:, 0]
Ya = analytic[:, 1]

fig = plt.figure(figsize=(24, 6))
fig.suptitle(
    "Time evolution of acceleration: comparison between analytical solution, explicit Euler, and RK2 schemes for differents time steps",
    fontsize=16,
)
gs = GridSpec(1, 3, figure=fig)

for i, (dt, tag) in enumerate(cases):
    # Load Euler
    euler = np.loadtxt(f"data/accel_euler_{tag}.dat", delimiter=";")
    Xe, Ye = euler[:, 0], euler[:, 1]

    # Load RK2
    rk2 = np.loadtxt(f"data/accel_rk2_{tag}.dat", delimiter=";")
    Xr, Yr = rk2[:, 0], rk2[:, 1]

    ax = fig.add_subplot(gs[0, i])
    ax.plot(Xe, Ye, lw=2, marker=".", label="Euler")
    ax.plot(Xr, Yr, lw=2, marker=".", label="RK2")
    ax.plot(Xa, Ya, lw=2, label="Analytic")

    ax.set_title(f"dt = {dt:.4f} s")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("v (m/s)")
    ax.grid()
    ax.legend()

plt.tight_layout()
plt.savefig("figures/Acceleration_schemes_comparison.png")
plt.show()


# %%
Error = np.loadtxt("data/Error_euler.dat", delimiter=";")
Xe = Error[:, 0]
Ye = Error[:, 1]

slope, intercept = np.polyfit(np.log(Xe), np.log(Ye), 1)
print("Order of convergence (Euler): ", slope)

plt.figure(figsize=(8, 8))
plt.plot(Xe, Ye, lw=2, marker="+", label="Euler")
plt.xlabel("dt (s)")
plt.ylabel("Max velocity error (m/s)")
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title(f"Log–log convergence of the explicit Euler method (slope = {slope:.2f})")
plt.legend()
plt.draw()
plt.savefig("figures/Error_Euler.png")


# %%
Error = np.loadtxt("data/Error_rk2.dat", delimiter=";")
Xe = Error[:, 0]
Ye = Error[:, 1]

# order of convergence
slope, _ = np.polyfit(np.log(Xe), np.log(Ye), 1)
print("Order of convergence (RK2): ", slope)

plt.figure(figsize=(8, 8))
plt.plot(Xe, Ye, lw=2, marker="v", label="RK2")
plt.xlabel("dt (s)")
plt.ylabel("Max velocity error (m/s)")
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title(f"Log–log convergence of the RK2 method (slope = {slope:.2f})")
plt.legend()
plt.draw()
plt.savefig("figures/Error_RK2.png")
plt.show()



