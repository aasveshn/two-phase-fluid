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
FRAME_SKIP = 1  # каждый N-й шаг 

if not os.path.exists(GIF_DIR):
    os.makedirs(GIF_DIR)

def load_mesh():
    
    mesh_path = os.path.join(INPUT_DIR, "mesh.txt")
    cells_map = []
    with open(mesh_path, 'r') as f:
        line = f.readline().split()
        if not line: return None
        nx, ny, num_cells, _ = map(int, line)
        for _ in range(num_cells):
            l = f.readline().split()
            cells_map.append((int(l[1]), int(l[2])))
    return nx, ny, cells_map

def line_to_array(line):

    return np.fromstring(line, sep=' ')

def get_grid(data, cells_map, nx, ny):
    
    grid = np.full((ny, nx), np.nan)
    for idx, (i, j) in enumerate(cells_map):
        grid[j, i] = data[idx]
    return grid

def data_generator(filenames, var_type, cells_map, nx, ny):
   
    
    handles = [open(os.path.join(OUTPUT_DIR, f), 'r') for f in filenames]
    
    try:
        count = 0
        while True:
            
            lines = [h.readline() for h in handles]
            
           
            if not any(lines) or any(l.strip() == "" for l in lines):
                break
            
            
            if count % FRAME_SKIP != 0:
                count += 1
                continue
            
            
            arrays = [line_to_array(l) for l in lines]
            
           
            if var_type == "pressure":
                
                res = arrays[0] / 1e5 
            elif var_type == "alpha":
                
                res = arrays[0]
            elif var_type == "density":
                
                a1, ro1, ro2 = arrays
                res = a1 * ro1 + (1.0 - a1) * ro2
            elif var_type == "velocity":
                
                a1, ro1, u1, v1, ro2, u2, v2 = arrays
                m1 = a1 * ro1
                m2 = (1.0 - a1) * ro2
                u_mix = (m1 * u1 + m2 * u2) / (m1 + m2)
                v_mix = (m1 * v1 + m2 * v2) / (m1 + m2)
                res = np.sqrt(u_mix**2 + v_mix**2)
            
            yield get_grid(res, cells_map, nx, ny)
            count += 1
    finally:
        for h in handles:
            h.close()

def create_animation(var_type, title, filenames, vmin=None, vmax=None, is_log=False):
    print(f"Обработка: {title}...")
    nx, ny, cells_map = load_mesh()
    extent = [0, nx * DX, 0, ny * DY]
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    
    
    gen = data_generator(filenames, var_type, cells_map, nx, ny)
    first_frame = next(gen)
    
    if is_log:
        
        v_min = vmin if vmin else max(np.nanmin(first_frame), 1e-5)
        v_max = vmax if vmax else np.nanmax(first_frame)
        norm = colors.LogNorm(vmin=v_min, vmax=v_max)
    else:
        norm = colors.Normalize(vmin=vmin, vmax=vmax)

    im = ax.imshow(first_frame, origin='lower', cmap='jet', extent=extent, norm=norm)
    plt.colorbar(im, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("x, м")
    ax.set_ylabel("y, м")
    ax.set_aspect('equal')

    def update(frame_grid):
        im.set_array(frame_grid)
        return [im]

   
    ani = FuncAnimation(fig, update, frames=gen, cache_frame_data=False, save_count=1000)
    
    save_path = os.path.join(GIF_DIR, f"{var_type}.gif")
    ani.save(save_path, writer='pillow', fps=FPS,dpi=100)
    plt.close(fig)
    print(f"Сохранено: {save_path}\n")



if __name__ == "__main__":
    
    create_animation("pressure", "Давление, бар", ["P1.txt"], is_log=True)

    
    create_animation("density", r"Плотность смеси, кг/м$^3$", ["a1.txt", "ro1.txt", "ro2.txt"], is_log=True)

    
    create_animation("velocity", "Модуль скорости смеси, м/с", 
                     ["a1.txt", "ro1.txt", "u1.txt", "v1.txt", "ro2.txt", "u2.txt", "v2.txt"], 
                     vmin=0, vmax=150, is_log=False)


    create_animation("alpha", "Объемная доля газа", ["a1.txt"], vmin=0, vmax=1, is_log=False)

    print("Готово")