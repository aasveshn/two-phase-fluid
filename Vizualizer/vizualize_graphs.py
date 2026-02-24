## Построение графиков в разных окнах
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import os

INPUT_DIR = "../input"
OUTPUT_DIR = "../output"

DX = 0.0004
DY = 0.0005

def load_data(filename):
    path = os.path.join(OUTPUT_DIR, filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Файл не найден: {path}")
    return np.loadtxt(path)

def plot_physical_results_log():
    mesh_path = os.path.join(INPUT_DIR, "mesh.txt")
    if not os.path.exists(mesh_path):
        print(f"Файл сетки не найден: {mesh_path}")
        return

    cells_map = []
    with open(mesh_path, 'r') as f:
        header = f.readline().split()
        if not header: 
            return
        nx_total, ny_total, num_cells, _ = map(int, header)
        for _ in range(num_cells):
            line = f.readline().split()
            cells_map.append((int(line[1]), int(line[2])))

    try:
        a1  = load_data("a1.txt")
        p1  = load_data("P1.txt")
        ro1 = load_data("ro1.txt")
        u1  = load_data("u1.txt")
        v1  = load_data("v1.txt")
        ro2 = load_data("ro2.txt")
        u2  = load_data("u2.txt")
        v2  = load_data("v2.txt")
    except Exception as e:
        print(f"Ошибка загрузки данных: {e}")
        return

    a2 = 1.0 - a1
    m1 = a1 * ro1
    m2 = a2 * ro2
    m_sum = m1 + m2
    
    p_bar = p1 / 1e5
    rho_mix = m_sum
    
    u_mix = (m1 * u1 + m2 * u2) / m_sum
    v_mix = (m1 * v1 + m2 * v2) / m_sum
    vel_mag = np.sqrt(u_mix**2 + v_mix**2)

    def to_grid(data_array):
        grid = np.full((ny_total, nx_total), np.nan)
        for idx, (i, j) in enumerate(cells_map):
            grid[j, i] = data_array[idx]
        return grid

    grid_p = to_grid(p_bar)
    grid_v = to_grid(vel_mag)
    grid_r = to_grid(rho_mix)
    grid_a = to_grid(a1)

    physical_extent = [0, nx_total * DX, 0, ny_total * DY]


    plt.figure(figsize=(8, 6))

    p_min = np.nanmin(grid_p)
    p_max = np.nanmax(grid_p)
    if p_min <= 0:
        p_min = 1e-6

    norm_p = colors.LogNorm(vmin=p_min, vmax=p_max)
    im1 = plt.imshow(grid_p, origin='lower', cmap='jet',
                 extent=physical_extent, norm=norm_p)

    cbar1 = plt.colorbar(im1)


    log_min = int(np.floor(np.log10(p_min)))
    log_max = int(np.ceil(np.log10(p_max)))

    ticks = []
    for n in range(log_min, log_max + 1):
        for m in range(1, 10):
            value = m * 10**n
            if p_min <= value <= p_max:
                ticks.append(value)

    cbar1.set_ticks(ticks)
    cbar1.set_ticklabels([f"{t:.1f}".rstrip('0').rstrip('.') for t in ticks])

    plt.title("Давление , бар")
    plt.xlabel("x, м")
    plt.ylabel("y, м")
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()


    plt.figure(figsize=(8, 6))

    r_min = np.nanmin(grid_r)
    r_max = np.nanmax(grid_r)
    if r_min <= 0:
        r_min = 1e-6

    norm_r = colors.LogNorm(vmin=r_min, vmax=r_max)
    im2 = plt.imshow(grid_r, origin='lower', cmap='jet',
                    extent=physical_extent, norm=norm_r)

    cbar2 = plt.colorbar(im2)

    log_min = int(np.floor(np.log10(r_min)))
    log_max = int(np.ceil(np.log10(r_max)))

    ticks = []
    for n in range(log_min, log_max + 1):
        for m in range(1, 10):
            value = m * 10**n
            if r_min <= value <= r_max:
                ticks.append(value)

    cbar2.set_ticks(ticks)
    cbar2.set_ticklabels([f"{t:.1f}".rstrip('0').rstrip('.') for t in ticks])

    plt.title(r"Плотность смеси , кг/м$^3$")
    plt.xlabel("x, м")
    plt.ylabel("y, м")
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()


    plt.figure(figsize=(8, 6))
    im3 = plt.imshow(grid_v, origin='lower', cmap='jet',
                     extent=physical_extent, vmin=0, vmax=150)
    plt.title("Модуль скорости смеси, м/с")
    plt.xlabel("x, м")
    plt.ylabel("y, м")
    plt.gca().set_aspect('equal')
    plt.colorbar(im3)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 6))
    im4 = plt.imshow(grid_a, origin='lower', cmap='jet',
                     extent=physical_extent, vmin=0, vmax=1)
    plt.title(r"Объемная доля газа ")
    plt.xlabel("x, м")
    plt.ylabel("y, м")
    plt.gca().set_aspect('equal')
    plt.colorbar(im4)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_physical_results_log()