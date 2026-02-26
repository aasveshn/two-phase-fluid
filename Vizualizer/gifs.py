import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.animation import FuncAnimation
import os

INPUT_DIR = "../input"
OUTPUT_DIR = "../output"
GIF_DIR = "gifs"

DX, DY = 0.0004, 0.0005
FPS = 10
FRAME_SKIP = 1  

if not os.path.exists(GIF_DIR):
    os.makedirs(GIF_DIR)

def load_mesh():
    mesh_path = os.path.join(INPUT_DIR, "mesh.txt")
    if not os.path.exists(mesh_path): return None
    cells_map = []
    with open(mesh_path, 'r') as f:
        header = f.readline().split()
        if not header: return None
        nx, ny, num_cells, _ = map(int, header)
        for _ in range(num_cells):
            l = f.readline().split()
            cells_map.append((int(l[1]), int(l[2])))
    return nx, ny, cells_map

def get_grid(data, cells_map, nx, ny):
    grid = np.full((ny, nx), np.nan)
    coords = np.array(cells_map)
    grid[coords[:, 1], coords[:, 0]] = data
    return grid

def process_line_to_physics(lines, var_type):
    arrays = [np.fromstring(l, sep=' ') for l in lines]
    
    if len(arrays[0]) == 0: return None

    if var_type == "pressure":
        return arrays[0] / 1e5
    elif var_type == "alpha":
        return arrays[0]
    elif var_type == "density":
        a1, ro1, ro2 = arrays
        return a1 * ro1 + (1.0 - a1) * ro2
    elif var_type == "velocity":
        a1, ro1, u1, v1, ro2, u2, v2 = arrays
        m1, m2 = a1 * ro1, (1.0 - a1) * ro2
        m_sum = m1 + m2
        u_mix = np.divide(m1*u1 + m2*u2, m_sum, out=np.zeros_like(m_sum), where=m_sum>1e-12)
        v_mix = np.divide(m1*v1 + m2*v2, m_sum, out=np.zeros_like(m_sum), where=m_sum>1e-12)
        return np.sqrt(u_mix**2 + v_mix**2)
    return None

def create_animation(var_type, title, filenames, vmin=None, vmax=None, is_log=False):
    print(f"Генерация GIF: {var_type}...")
    nx, ny, cells_map = load_mesh()
    extent = [0, nx * DX, 0, ny * DY]

    handles = [open(os.path.join(OUTPUT_DIR, f), 'r') for f in filenames]
    
    frames_list = []
    
    for count, lines in enumerate(zip(*handles)):
        data = process_line_to_physics(lines, var_type)
        
        if data is not None:
            if count % FRAME_SKIP == 0:
                frames_list.append(data)
            last_step_data = data 

    if count % FRAME_SKIP != 0:
        frames_list.append(last_step_data)

    for h in handles: h.close()

    if not frames_list:
        print(f"Данные для {var_type} не найдены.")
        return

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    
    all_vals = np.concatenate(frames_list)
    if is_log:
        v_min = vmin if vmin else max(np.nanmin(all_vals), 1e-5)
        v_max = vmax if vmax else np.nanmax(all_vals)
        norm = colors.LogNorm(vmin=v_min, vmax=v_max)
    else:
        v_min = vmin if vmin is not None else np.nanmin(all_vals)
        v_max = vmax if vmax is not None else np.nanmax(all_vals)
        norm = colors.Normalize(vmin=v_min, vmax=v_max)

    curr_grid = get_grid(frames_list[0], cells_map, nx, ny)
    im = ax.imshow(curr_grid, origin='lower', cmap='jet', extent=extent, norm=norm)
    plt.colorbar(im, ax=ax)
    ax.set_title(title)
    ax.set_aspect('equal')

    def update(frame_data):
        grid = get_grid(frame_data, cells_map, nx, ny)
        im.set_array(grid)
        return [im]

    ani = FuncAnimation(fig, update, frames=frames_list, blit=True, cache_frame_data=False)
    
    save_path = os.path.join(GIF_DIR, f"{var_type}.gif")
    ani.save(save_path, writer='pillow', fps=FPS)
    plt.close(fig)
    print(f"Выполненно\n")

if __name__ == "__main__":
    create_animation("pressure", "Давление, бар", ["P1.txt"], is_log=True)
    create_animation("density", "Плотность смеси, кг/м3", ["a1.txt", "ro1.txt", "ro2.txt"], is_log=True)
    create_animation("velocity", "Модуль скорости смеси, м/с", 
                     ["a1.txt", "ro1.txt", "u1.txt", "v1.txt", "ro2.txt", "u2.txt", "v2.txt"], 
                     vmin=0, vmax=150, is_log=False)
    create_animation("alpha", "Объемная доля газа", ["a1.txt"], vmin=0, vmax=1, is_log=False)