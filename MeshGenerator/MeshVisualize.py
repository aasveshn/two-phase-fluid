import matplotlib.pyplot as plt
import matplotlib.patches as patches

def visualize_mesh(filename, dx=0.01, dy=0.01):
    cells = []
    faces = []

    with open(filename, 'r') as f:
        header = f.readline().split()
        if not header: return
        nx_total, ny_total, num_cells, num_faces = map(int, header)

        for _ in range(num_cells):
            line = f.readline().split()
            # ID, i, j, f_left, f_right, f_bottom, f_top
            cells.append({
                'id': int(line[0]),
                'i': int(line[1]),
                'j': int(line[2])
            })

        for _ in range(num_faces):
            line = f.readline().split()
            # ID, L_Cell, R_Cell, isVertical, Type
            faces.append({
                'id': int(line[0]),
                'L': int(line[1]),
                'R': int(line[2]),
                'isVert': int(line[3]),
                'type': int(line[4])
            })

    fig, ax = plt.subplots(figsize=(12, 6))
    
    for cell in cells:
        x = cell['i'] * dx
        y = cell['j'] * dy
        rect = patches.Rectangle((x, y), dx, dy, linewidth=0, edgecolor='none', facecolor='skyblue', alpha=0.3)
        ax.add_patch(rect)

    colors = {
        0: 'lightgray', # внутренняя 
        1: 'black',     # стенка 
        2: 'blue',      # вход 
        3: 'red'        # выход 
    }
    
    
    for f in faces:
        
        if f['isVert']:
            target_cid = f['R'] if f['R'] != -1 else f['L']
            target_cell = next(c for c in cells if c['id'] == target_cid)
            
            idx_i = target_cell['i'] if f['R'] == target_cid else target_cell['i'] + 1
            x = [idx_i * dx, idx_i * dx]
            y = [target_cell['j'] * dy, (target_cell['j'] + 1) * dy]
        
        else:
            target_cid = f['R'] if f['R'] != -1 else f['L']
            target_cell = next(c for c in cells if c['id'] == target_cid)
            
            idx_j = target_cell['j'] if f['R'] == target_cid else target_cell['j'] + 1
            x = [target_cell['i'] * dx, (target_cell['i'] + 1) * dx]
            y = [idx_j * dy, idx_j * dy]

        lw = 2 if f['type'] != 0 else 0.5 
        ax.plot(x, y, color=colors[f['type']], linewidth=lw)


    ax.set_aspect('equal')
    ax.set_title(f'Визуализация сетки: {num_cells} ячеек, {num_faces} граней')
    ax.set_xlabel('X, м')
    ax.set_ylabel('Y, м')
    
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='black', lw=2),
                    Line2D([0], [0], color='blue', lw=2),
                    Line2D([0], [0], color='red', lw=2),
                    Line2D([0], [0], color='lightgray', lw=1)]
    ax.legend(custom_lines, ['Стенка', 'Вход', 'Выход', 'Внутренняя'])

    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

if __name__ == "__main__":
    visualize_mesh("mesh.txt", dx=0.01, dy=0.01)
