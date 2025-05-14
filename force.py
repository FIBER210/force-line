import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import numpy as np
import re

def point_to_segment_distance(px, py, x1, y1, x2, y2):
    """
    Calculate the minimum distance from point (px, py) to the line segment defined by (x1, y1) and (x2, y2).
    """
    # Vector from point 1 to point 2
    dx = x2 - x1
    dy = y2 - y1

    if dx == 0 and dy == 0:
        # The segment is a single point
        return np.hypot(px - x1, py - y1)

    # Parameter t of the projection of point p onto the line defined by the segment
    t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy)

    if t < 0:
        # Closest to point 1
        closest_x, closest_y = x1, y1
    elif t > 1:
        # Closest to point 2
        closest_x, closest_y = x2, y2
    else:
        # Projection falls on the segment
        closest_x = x1 + t * dx
        closest_y = y1 + t * dy

    return np.hypot(px - closest_x, py - closest_y)

def point_to_polyline_distance(px, py, polyline_points):
    """
    Calculate the minimum distance from point (px, py) to a polyline defined by a list of points [(x0, y0), (x1, y1), ...].
    """
    min_dist = float('inf')
    for i in range(len(polyline_points) - 1):
        x1, y1 = polyline_points[i]
        x2, y2 = polyline_points[i + 1]
        dist = point_to_segment_distance(px, py, x1, y1, x2, y2)
        if dist < min_dist:
            min_dist = dist
    return min_dist

def segments_intersect(p1, p2, q1, q2):
    """
    Check if line segment p1-p2 intersects with line segment q1-q2.
    p1, p2, q1, q2 are tuples (x, y).
    """
    def orientation(a, b, c):
        # Returns the orientation of the triplet (a, b, c)
        # 0 -> colinear, 1 -> clockwise, 2 -> counterclockwise
        val = (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1])
        if abs(val) < 1e-10:
            return 0
        return 1 if val > 0 else 2

    def on_segment(a, b, c):
        # Check if point b lies on segment a-c
        if min(a[0], c[0]) - 1e-10 <= b[0] <= max(a[0], c[0]) + 1e-10 and \
           min(a[1], c[1]) - 1e-10 <= b[1] <= max(a[1], c[1]) + 1e-10:
            return True
        return False

    o1 = orientation(p1, p2, q1)
    o2 = orientation(p1, p2, q2)
    o3 = orientation(q1, q2, p1)
    o4 = orientation(q1, q2, p2)

    # General case
    if o1 != o2 and o3 != o4:
        return True

    # Special Cases
    # p1, p2 and q1 are colinear and q1 lies on segment p1p2
    if o1 == 0 and on_segment(p1, q1, p2):
        return True

    # p1, p2 and q2 are colinear and q2 lies on segment p1p2
    if o2 == 0 and on_segment(p1, q2, p2):
        return True

    # q1, q2 and p1 are colinear and p1 lies on segment q1q2
    if o3 == 0 and on_segment(q1, p1, q2):
        return True

    # q1, q2 and p2 are colinear and p2 lies on segment q1q2
    if o4 == 0 and on_segment(q1, p2, q2):
        return True

    return False

def segment_intersects_polyline(seg_start, seg_end, polyline_points):
    """
    Check if the line segment defined by seg_start and seg_end intersects with the polyline.
    seg_start, seg_end: tuples (x, y)
    polyline_points: list of tuples [(x0, y0), (x1, y1), ...]
    Returns True if any segment of the polyline intersects with the given segment.
    """
    for i in range(len(polyline_points) - 1):
        poly_start = polyline_points[i]
        poly_end = polyline_points[i + 1]
        if segments_intersect(seg_start, seg_end, poly_start, poly_end):
            return True
    return False

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


def extract_multiple_polylines(file_path, skip_header=0):
    polylines = {}  # 存储 {编号: [(x1, y1), (x2, y2), ...]}
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_number, line in enumerate(f, start=1):
                if line_number <= skip_header:
                    continue  # Skip header lines
                # Split by comma (normal or full-width) or whitespace
                values = re.split(r'[,\s，]+', line.strip())
                if len(values) >= 3:  # 至少 x, y, 编号
                    try:
                        x = float(values[0])
                        y = float(values[1])
                        polyline_id = int(values[2])  # 假设编号是整数
                        if polyline_id not in polylines:
                            polylines[polyline_id] = []
                        polylines[polyline_id].append((x, y))
                    except ValueError:
                        print(f"Warning: Line {line_number} has invalid data and will be skipped.")
                else:
                    print(f"Warning: Line {line_number} does not have enough data and will be skipped.")
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")
    return polylines

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
def get_next_point(current_point, step_size, direction=1):
    sx, sy = current_point.get_stress_vector()
    magnitude = np.hypot(sx, sy)

    unit_sx = sx / magnitude
    unit_sy = sy / magnitude


    # 将方向向量转换为复数并乘以 direction
    vec = complex(unit_sx, unit_sy) * direction

    next_x = current_point.x + vec.real * step_size
    next_y = current_point.y + vec.imag * step_size
    return Point(next_x, next_y)

def trace_force_path(start_point, points_collection, delaunay, coords, step_size, max_steps,edge):
    for simplex in delaunay.simplices:
        triangle_points = coords[simplex]
        if is_point_in_triangle(start_point, triangle_points):
            get_stress_at_triangle(start_point, triangle_points, points_collection)
            break  # # 绘制起始点
    current_point=start_point
    path = [current_point]
    for step in range(max_steps):
        # 1. 计算下一个点（沿当前应力方向移动）
        next_point = get_next_point(current_point, step_size,1)
        print(f"Forward step {step+1}: Point({next_point.x:.2f}, {next_point.y:.2f})")

        if any(
                segment_intersects_polyline(
                    (current_point.x, current_point.y),
                    (next_point.x, next_point.y),
                    polyline
                )
                for polyline in edge.values()  # 遍历所有多段线
        ):
            print(f"路径在 ({next_point.x:.2f}, {next_point.y:.2f}) 与边界相交，终止计算")
            break
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
        next_point = get_next_point(current_point, step_size,-1)

        print(f"Forward step {step+1}: Point({next_point.x:.2f}, {next_point.y:.2f})")

        if any(
                segment_intersects_polyline(
                    (current_point.x, current_point.y),
                    (next_point.x, next_point.y),
                    polyline
                )
                for polyline in edge.values()  # 遍历所有多段线
        ):
            print(f"路径在 ({next_point.x:.2f}, {next_point.y:.2f}) 与边界相交，终止计算")
            break
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
    file_path = ("slice1.txt")
    edge = extract_multiple_polylines(file_path)

    # 提取坐标
    coords = np.array([(point.x, point.y) for point in points_collection.get_points()])

    # 使用 Delaunay 三角剖分生成三角形
    delaunay = Delaunay(coords)

    # 设置指定的起始点 (手动指定起始点坐标)
    start_points = [
        # Point(x=40, y=22), Point(x=40, y=21), Point(x=40, y=20),
        # Point(x=40, y=19), Point(x=40, y=18), Point(x=40, y=17),
        # Point(x=40, y=16), Point(x=40, y=15), Point(x=40, y=14),
        # Point(x=40, y=13), Point(x=40, y=12), Point(x=40, y=11),
        # Point(x=40, y=10), Point(x=40, y=9), Point(x=40, y=8),
        # Point(x=40, y=7), Point(x=40, y=6), Point(x=40, y=5),
        # Point(x=40, y=4), Point(x=40, y=3), Point(x=40, y=2),
        # Point(x=40, y=1), Point(x=40, y=0), Point(x=40, y=-1),
        # Point(x=40, y=-2), Point(x=40, y=-3), Point(x=40, y=-4),
        # Point(x=40, y=-5), Point(x=40, y=-6), Point(x=40, y=-7),
        # Point(x=40, y=-8), Point(x=40, y=-9), Point(x=40, y=-10),
        # Point(x=40, y=-11), Point(x=40, y=-12), Point(x=40, y=-13),
        # Point(x=40, y=-14), Point(x=40, y=-15), Point(x=40, y=-16),
        # Point(x=40, y=-17), Point(x=40, y=-18), Point(x=40, y=-19),
        # Point(x=40, y=-20), Point(x=40, y=-21), Point(x=40, y=-22)
        Point(40, 0),  # 起始点2
        # Point(0, 19),  # 起始点2
        # Point(0, 16), # 起始点3
        # Point(0, 17),  # 起始点2
        # Point(0, 20),  # 起始点3
        # Point(0, 21),  # 起始点2
        # Point(0, 22) , # 起始点3
        # Point(0, -22),  # 起始点3
        # Point(0, -21),  # 起始点3
        # Point(0, -16),  # 起始点3
        # Point(0, -18),  # 起始点3
        # Point(0, -17),  # 起始点3
        # Point(0, -20),  # 起始点3
        # Point(0, -19),  # 起始点3

    ]
    # start_point = Point(4, 5, 0, 0)  # 假设应力方向为 (0, 0)

    # 进行力流路径的迭代
    plt.figure(figsize=(10, 8))
    ax = plt.gca()

    step_size = 0.1  # 每次步进的距离
    max_steps = 500

    for i, start_point in enumerate(start_points):



        path=trace_force_path(start_point, points_collection, delaunay, coords, step_size, max_steps,edge)


        if len(path) > 1:
            # 路径线（红色）
            ax.plot([p.x for p in path], [p.y for p in path], 'r-', lw=2, label='应力流径')



    # 1. 绘制三角网格背景
    for simplex in delaunay.simplices:
        simplex = np.append(simplex, simplex[0])  # 使得三角形闭合
        ax.plot(coords[simplex][:, 0], coords[simplex][:, 1], 'gray', alpha=0.2)

    # 2. 绘制原始数据点
    ax.scatter(coords[:, 0], coords[:, 1], color='blue', s=10, alpha=0.5)

    for polyline_id, points in edge.items():
        x_coords, y_coords = zip(*points)  # 解压坐标点
        ax.plot(
            x_coords,
            y_coords,
            marker='o',  # 点标记样式
            linestyle='-',  # 线型
            color='blue',  # 颜色（可以改成不同颜色区分不同多段线）
            label=f'Path {polyline_id}'  # 添加图例标签
        )

    plt.show()



if __name__ == "__main__":
    main()
