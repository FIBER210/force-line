import trimesh
import numpy as np

def load_stl(file_path):
    # 使用 trimesh 加载 STL 文件
    return trimesh.load_mesh(file_path)

def get_bounding_box(mesh):
    # 获取网格的坐标范围
    min_coords = mesh.bounds[0]
    max_coords = mesh.bounds[1]
    return min_coords, max_coords

def get_cross_section(mesh, z_height):
    # 在指定的 Z 高度处获取截面
    plane_normal = [0, 0, 1]  # Z 平面的法向量
    plane_origin = [0, 0, z_height]  # Z 高度

    # 使用 trimesh 进行截面提取
    section = mesh.section(plane_origin=plane_origin, plane_normal=plane_normal)

    if section is None:
        return None

    # 获取截面轮廓
    return section


def save_to_txt(cross_section, output_file):
    # 投影为2D平面路径
    cross_section_2D, _ = cross_section.to_planar()

    vertices = np.array(cross_section_2D.vertices)
    entities = cross_section_2D.entities

    with open(output_file, 'w') as f:
        for line_index, entity in enumerate(entities):
            # 每个 entity 是一条线段，有两个端点索引
            if hasattr(entity, 'points'):
                # 对于 Line 类型，可以直接用点
                points = entity.discrete(vertices)
                for point in points:
                    f.write(f"{point[0]:.6f}, {point[1]:.6f}, {line_index}\n")


def main(stl_file, z_height, txt_file):
    # 加载STL文件
    mesh = load_stl(stl_file)

    # 获取坐标范围
    min_coords, max_coords = get_bounding_box(mesh)
    print(f"Bounding Box: Min {min_coords}, Max {max_coords}")

    # 获取指定z高度的截面
    cross_section = get_cross_section(mesh, z_height)
    if cross_section is not None:
        print(f"Cross section at z={z_height}: {cross_section}")
        # 保存为TXT文件
        save_to_txt(cross_section, txt_file)
        print(f"TXT file saved as: {txt_file}")
    else:
        print(f"No valid cross section found at z={z_height}")

# 示例调用
stl_file = 'hole.STL'
z_height = 0.5  # 指定的Z高度,单位是mm
txt_file = 'slice.txt'
main(stl_file, z_height, txt_file)
