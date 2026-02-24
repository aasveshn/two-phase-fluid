import numpy as np

# =================================================================
# =========== ГЕОМЕТРИЧЕСКИЕ ПАРАМЕТРЫ И РАЗРЕШЕНИЕ============
# =================================================================
dx = 0.0004 # шаг сетки по X
dy = 0.0005  # шаг сетки по Y

# параметры трубы
pipe_length = 0.150  # длина трубы
pipe_height = 0.005  # высота трубы

# параметры камеры
chamber_length = 0.050  # длина камеры
chamber_height = 0.02  # высота камеры

# типы граней 
TYPE_INTERNAL = 0 #внутреняя
TYPE_WALL     = 1 # стена (в трубе и слева в камере)
TYPE_INLET    = 2  # слева в трубе
TYPE_OUTLET   = 3  # справа в камере


nx_pipe = int(round(pipe_length / dx))
ny_pipe = int(round(pipe_height / dy))

nx_chamber = int(round(chamber_length / dx))
ny_chamber = int(round(chamber_height / dy))


nx_total = nx_pipe + nx_chamber
ny_total = max(ny_pipe, ny_chamber)


pipe_y_start = (ny_total - ny_pipe) // 2
pipe_y_end = pipe_y_start + ny_pipe


cell_map = np.full((nx_total, ny_total), -1, dtype=int)
cell_list = []

cell_id_counter = 0
for i in range(nx_total):
    for j in range(ny_total):
        is_fluid = False
        

        if i < nx_pipe:
            if pipe_y_start <= j < pipe_y_end:
                is_fluid = True
    
        else:
            if 0 <= j < ny_chamber:
                is_fluid = True
        
        if is_fluid:
            cell_map[i, j] = cell_id_counter
            cell_list.append((cell_id_counter, i, j))
            cell_id_counter += 1

num_cells = cell_id_counter


faces_list = [] # (id, L_cell, R_cell, is_vertical, type)
face_id_counter = 0


cell_faces = {cid: [-1, -1, -1, -1] for cid in range(num_cells)}


for i in range(nx_total + 1):
    for j in range(ny_total):
        L_cid = cell_map[i-1, j] if i > 0 else -1
        R_cid = cell_map[i, j] if i < nx_total else -1
        
        if L_cid == -1 and R_cid == -1: continue 
        
        f_type = TYPE_INTERNAL
        if L_cid == -1 or R_cid == -1:
            if i == 0: f_type = TYPE_INLET
            elif i == nx_total: f_type = TYPE_OUTLET
            else: f_type = TYPE_WALL
        
        fid = face_id_counter
        faces_list.append((fid, L_cid, R_cid, 1, f_type)) 
        
        if L_cid != -1: cell_faces[L_cid][1] = fid
        if R_cid != -1: cell_faces[R_cid][0] = fid
        face_id_counter += 1


for i in range(nx_total):
    for j in range(ny_total + 1):
        B_cid = cell_map[i, j-1] if j > 0 else -1 
        T_cid = cell_map[i, j] if j < ny_total else -1 
        
        if B_cid == -1 and T_cid == -1: continue
        
        f_type = TYPE_INTERNAL
        if B_cid == -1 or T_cid == -1:
            f_type = TYPE_WALL 
            
        fid = face_id_counter
        faces_list.append((fid, B_cid, T_cid, 0, f_type)) 
        
        if B_cid != -1: cell_faces[B_cid][3] = fid 
        if T_cid != -1: cell_faces[T_cid][2] = fid 
        face_id_counter += 1

num_faces = face_id_counter


with open("../mesh.txt", "w") as f:

    f.write(f"{nx_total} {ny_total} {num_cells} {num_faces}\n")
    
    
   
    for cid, i, j in cell_list:
        fids = cell_faces[cid]
        f.write(f"{cid} {i} {j} {fids[0]} {fids[1]} {fids[2]} {fids[3]}\n")
    
    for fid, L, R, is_vert, f_type in faces_list:
        f.write(f"{fid} {L} {R} {is_vert} {f_type}\n")

print(f"Сетка успешно сгенерирована!")
print(f"Ячеек: {num_cells}, Граней: {num_faces}")
print(f"Размер области: {nx_total} x {ny_total}")
