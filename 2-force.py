import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from sklearn.neighbors import NearestNeighbors

class Point:
    def __init__(self, x, y, stress_x=0.0, stress_y=0.0):
        self.x = x  # x 坐标
        self.y = y  # y 坐标
        self.stress_x = stress_x  # 第一主应力的 x 分量
        self.stress_y = stress_y  # 第一主应力的 y 分量

    def get_stress_vector(self):
        # 返回应力矢量（x 和 y 分量）
        return (self.stress_x, self.stress_y)

    def __str__(self):
        # 返回点的描述信息
        return f"Point(x={self.x}, y={self.y}, stress_vector=({self.stress_x}, {self.stress_y}))"

    def add_stress(self, add_stress_x, add_stress_y):

            self.stress_x += add_stress_x
            self.stress_y += add_stress_y

def find_point_by_xy(points, x, y):
    for point in points:
        if point.x == x and point.y == y:
            return point

class PointCollection:
    def __init__(self):
        self.points = {}  # 用字典存储，键为 (x, y)，值为 Point 对象

    def add_point(self, point):
        self.points[(point.x, point.y)] = point

    def get_point(self, x, y, tol=1e-6):
        for (px, py), p in self.points.items():
            if abs(px - x) < tol and abs(py - y) < tol:
                return p
        return None

    def get_points(self):
        return list(self.points.values())

    def to_array(self):
        return np.array([
            [p.x, p.y, p.stress_x, p.stress_y]
            for p in self.points.values()
        ])
# 解析 TXT 文件并构建点对象
def load_points_from_file(file_path):
    points = PointCollection()  # 使用集合存储点对象
    try:
        with open(file_path, 'r') as f:
            for line in f:
                values = line.strip().split()  # 假设每个值用空格分隔
                if len(values) >= 5:  # 确保每行至少有5个值
                    x, y, z, stress_x, stress_y = map(float, values[:5])  # 提取 x, y, 第一主应力
                    point = Point(x, y, stress_x, stress_y)
                    points.add_point(point)
    except FileNotFoundError:
        print(f"错误: 文件 {file_path} 未找到。")
    except ValueError:
        print("错误: 文件中的数据格式不正确，无法转换为浮点数。")
    except Exception as e:
        print(f"发生未知错误: {e}")
    return points

# 从给定的三角形网格中判断点是否在其中
def is_point_in_triangle(point, triangle_points, tol=1e-6):
    x, y = point.x, point.y
    x1, y1 = triangle_points[0]
    x2, y2 = triangle_points[1]
    x3, y3 = triangle_points[2]

    area_orig = abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    area1 = abs(x * (y2 - y3) + x2 * (y3 - y) + x3 * (y - y2))
    area2 = abs(x1 * (y - y3) + x * (y3 - y1) + x3 * (y1 - y))
    area3 = abs(x1 * (y2 - y) + x2 * (y - y1) + x * (y1 - y2))

    return abs(area_orig - (area1 + area2 + area3)) < tol


# 获取三角形中的点的应力值
def get_stress_at_triangle(point, triangle_points, points_collection):
    """
    如果 point 在 triangle_points 所定义的三角形内：
    - 使用三角形三个点的应力值进行插值
    - 并将插值结果赋值给 point 的 stress_x 和 stress_y

    参数:
        point (Point): 目标点
        triangle_points (list of tuple): 三角形顶点 [(x1, y1), (x2, y2), (x3, y3)]
        points_collection (PointCollection): 所有点集合

    返回:
        bool: 如果插值成功返回 True，否则返回 False
    """
    if not is_point_in_triangle(point, triangle_points):
        return False

    # 获取三角形三个顶点
    (x1, y1), (x2, y2), (x3, y3) = triangle_points

    # 查找三个顶点对应的点对象
    p1 = find_point_by_xy(points_collection.get_points(), x1, y1)
    p2 = find_point_by_xy(points_collection.get_points(), x2, y2)
    p3 = find_point_by_xy(points_collection.get_points(), x3, y3)

    if not p1 or not p2 or not p3:
        print("警告：三角形顶点中有点未找到，应力插值失败。")
        return False

    # 使用重心坐标计算插值权重
    x, y = point.x, point.y
    detT = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
    if abs(detT) < 1e-10:
        print("警告：三角形退化，面积为0。")
        return False

    # 重心系数
    lambda1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / detT
    lambda2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / detT
    lambda3 = 1 - lambda1 - lambda2

    # 插值应力
    stress_x = lambda1 * p1.stress_x + lambda2 * p2.stress_x + lambda3 * p3.stress_x
    stress_y = lambda1 * p1.stress_y + lambda2 * p2.stress_y + lambda3 * p3.stress_y

    # 赋值给目标点
    point.stress_x = stress_x
    point.stress_y = stress_y

    return True


# 计算下一点的坐标
def get_next_point(current_point, step_size):
    sx, sy = current_point.get_stress_vector()
    next_x = current_point.x + sx * step_size
    next_y = current_point.y + sy * step_size
    return Point(next_x, next_y)

def trace_force_path(start_point, points_collection, delaunay, coords, step_size, max_steps, ):
    current_point=start_point
    path = [current_point]
    for step in range(max_steps):
        # 1. 计算下一个点（沿当前应力方向移动）
        next_x = current_point.x + current_point.stress_x * step_size
        next_y = current_point.y + current_point.stress_y * step_size
        next_point = Point(next_x, next_y)  # 新点初始应力为0

        # 2. 判断新点是否在任何三角形内
        point_in_mesh = False
        for simplex in delaunay.simplices:
            triangle = coords[simplex]  # 获取三角形顶点坐标

            if is_point_in_triangle(next_point, triangle):
                # 3. 如果在三角形内 -> 计算插值应力
                if get_stress_at_triangle(next_point, triangle, points_collection):
                    point_in_mesh = True
                    break  # 找到所属三角形后立即跳出循环

        # 4. 根据判断结果决定是否继续
        if not point_in_mesh:
            print(f"路径在 ({next_point.x:.2f}, {next_point.y:.2f}) 离开网格，终止计算")
            break

        # 5. 将有效点加入路径并更新当前点
        path.append(next_point)
        current_point = next_point

        # 更新路径数据
        path_x = [p.x for p in path]
        path_y = [p.y for p in path]


    current_point = start_point
    for step in range(max_steps):
        # 1. 计算下一个点（沿当前应力方向移动）
        next_x = current_point.x - current_point.stress_x * step_size
        next_y = current_point.y - current_point.stress_y * step_size
        next_point = Point(next_x, next_y)  # 新点初始应力为0

        # 2. 判断新点是否在任何三角形内
        point_in_mesh = False
        for simplex in delaunay.simplices:
            triangle = coords[simplex]  # 获取三角形顶点坐标

            if is_point_in_triangle(next_point, triangle):
                # 3. 如果在三角形内 -> 计算插值应力
                if get_stress_at_triangle(next_point, triangle, points_collection):
                    point_in_mesh = True
                    break  # 找到所属三角形后立即跳出循环

        # 4. 根据判断结果决定是否继续
        if not point_in_mesh:
            print(f"路径在 ({next_point.x:.2f}, {next_point.y:.2f}) 离开网格，终止计算")
            break

        # 5. 将有效点加入路径并更新当前点
        path.insert(0, next_point)
        current_point = next_point

        # 更新路径数据
        path_x = [p.x for p in path]
        path_y = [p.y for p in path]

    return path
# 主程序
def main():
    file_path = ("4-28foce_processed.txt")  # 假设 txt 文件名为 points_data.txt
    points_collection = load_points_from_file(file_path)
    data = points_collection.to_array()
    positions =  data[:, :2]  # (x, y)
    directions =  data[:, 2:]  # (σx, σy)
    n_points = len(positions)

    # ⚙️ 构建稀疏矩阵系统 A * s = b
    k = 5  # 每个点找 k 个邻居
    neighbors = NearestNeighbors(n_neighbors=k).fit(positions)
    _, indices = neighbors.kneighbors(positions)

    A = lil_matrix((n_points, n_points))
    b = np.zeros(n_points)


    A[0, :] = 0
    A[0, 0] = 1
    b[0] = 0

    for i in range(1, n_points):  # 从第2个点开始
        xi = positions[i]
        ui = directions[i]
        norm = np.linalg.norm(ui)
        if norm < 1e-6:
            A[i, i] = 1
            b[i] = 0
            continue
        ui /= norm
        A[i, i] = k
        for j in indices[i][1:]:
            xj = positions[j]
            dx = xj - xi
            proj = np.dot(dx, ui)
            A[i, j] = -1
            b[i] += proj

    # 求解线性系统
    from scipy.sparse import csr_matrix
    s = spsolve(csr_matrix(A), b)

    # 可视化标量场（散点图）
    plt.figure(figsize=(6, 5))
    plt.scatter(positions[:, 0], positions[:, 1], c=s, cmap='viridis', s=40)
    plt.colorbar(label='Scalar Field s(x)')
    plt.title('Scalar Field Visualization (scatter)')
    plt.axis('equal')
    plt.show()

    # 可视化插值（等高线图）
    x_grid = np.linspace(positions[:, 0].min(), positions[:, 0].max(), 100)
    y_grid = np.linspace(positions[:, 1].min(), positions[:, 1].max(), 100)
    grid_x, grid_y = np.meshgrid(x_grid, y_grid)
    grid_s = griddata(positions, s, (grid_x, grid_y), method='cubic')

    plt.figure(figsize=(6, 5))
    cs = plt.contourf(grid_x, grid_y, grid_s, levels=100, cmap='viridis')
    plt.colorbar(cs, label='s(x)')
    plt.title('Interpolated Scalar Field (contour)')
    plt.axis('equal')
    plt.show()

    # 可选：叠加梯度箭头（方向场）
    gx, gy = np.gradient(grid_s)
    plt.figure(figsize=(6, 5))
    plt.contourf(grid_x, grid_y, grid_s, levels=50, cmap='viridis')
    plt.quiver(grid_x, grid_y, gx, gy, color='white', scale=50)
    plt.title('Gradient Field ∇s(x)')
    plt.axis('equal')
    plt.show()

    # # 提取坐标
    # coords = np.array([(point.x, point.y) for point in points_collection.get_points()])
    #
    # # 使用 Delaunay 三角剖分生成三角形
    # delaunay = Delaunay(coords)
    #
    # # 设置指定的起始点 (手动指定起始点坐标)
    # start_points = [
    #     Point(45, 20.0),  # 起始点1
    #     # Point(45, 15),  # 起始点2
    #     # Point(45, 10), # 起始点3
    #     # Point(45, 5.0) , # 起始点3
    #     # Point(45, 0.05),  # 起始点3
    #     # Point(45, -5.0),  # 起始点3
    #     # Point(45, -10.0),  # 起始点3
    #     # Point(45, -15.0),  # 起始点3
    #     # Point(45, -20.0),  # 起始点3
    #
    # ]
    # # start_point = Point(4, 5, 0, 0)  # 假设应力方向为 (0, 0)
    #
    # # 进行力流路径的迭代
    # plt.figure(figsize=(10, 8))
    # ax = plt.gca()
    #
    # step_size = 0.1  # 每次步进的距离
    # max_steps = 10000
    #
    # for i, start_point in enumerate(start_points):
    #
    #     for simplex in delaunay.simplices:
    #         triangle_points = coords[simplex]
    #         if is_point_in_triangle(start_point, triangle_points):
    #             get_stress_at_triangle(start_point, triangle_points, points_collection)
    #             break    # # 绘制起始点
    #
    #     path=trace_force_path(start_point, points_collection, delaunay, coords, step_size, max_steps)
    #
    #
    #     if len(path) > 1:
    #         # 路径线（红色）
    #         ax.plot([p.x for p in path], [p.y for p in path], 'r-', lw=2, label='应力流径')
    #
    #
    #
    # # 1. 绘制三角网格背景
    # for simplex in delaunay.simplices:
    #     simplex = np.append(simplex, simplex[0])  # 使得三角形闭合
    #     ax.plot(coords[simplex][:, 0], coords[simplex][:, 1], 'gray', alpha=0.2)
    #
    # # 2. 绘制原始数据点
    # ax.scatter(coords[:, 0], coords[:, 1], color='blue', s=10, alpha=0.5)
    #
    # # 3. 绘制完整路径
    #
    #
    #     # 路径点颜色映射应力大小
    #     # stresses = [np.sqrt(p.stress_x ** 2 + p.stress_y ** 2) for p in path]
    #     # sc = ax.scatter(
    #     #     [p.x for p in path], [p.y for p in path],
    #     #     c=stresses, cmap='coolwarm', s=50,
    #     #     edgecolors='k', linewidths=0.5,
    #     #     label='应力值'
    #     # )
    #     # plt.colorbar(sc, label='应力大小 (MPa)')
    #
    #     # # 标记关键点
    #     # ax.scatter(path[0].x, path[0].y, color='lime', s=150, marker='*', label='起始点')
    #     # ax.scatter(path[-1].x, path[-1].y, color='red', s=150, marker='X', label='终止点')
    #
    #     # 添加方向箭头（稀疏显示）
    #     # for i in range(0, len(path), len(path) // 10 + 1):  # 约10个箭头
    #     #     p = path[i]
    #     #     ax.quiver(
    #     #         p.x, p.y, p.stress_x, p.stress_y,
    #     #         color='black', scale=30, width=0.004,
    #     #         headwidth=4, alpha=0.7
    #     #     )
    # #     # 4. 图形标注
    # # ax.set_xlabel('X坐标 (mm)', fontsize=12)
    # # ax.set_ylabel('Y坐标 (mm)', fontsize=12)
    # # ax.set_title('应力流径分析结果', fontsize=14)
    # # ax.legend(loc='upper right')
    # # ax.grid(True, linestyle='--', alpha=0.3)
    # # plt.tight_layout()
    # plt.show()



if __name__ == "__main__":
    main()
