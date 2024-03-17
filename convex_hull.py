import tkinter as tk
from tkinter import Canvas, ttk
from functools import cmp_to_key
import numpy as np
from scipy.spatial import ConvexHull
import time

class Point:
    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y

def direction(p, q, r):
    return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)

def distance_sq(p1, p2):
    return (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2

def jarvis_march(points):
    a = min(points, key=lambda point: point.x)
    index = points.index(a)

    l = index
    result = []
    result.append(a)
    while True:
        q = (l + 1) % len(points)
        for i in range(len(points)):
            if i == l:
                continue
            d = direction(points[l], points[i], points[q])
            if d > 0 or (d == 0 and distance_sq(points[i], points[l]) > distance_sq(points[q], points[l])):
                q = i
        l = q
        if l == index:
            break
        result.append(points[q])

    return result

def monotone_chain(points):
    def orientation(p, q, r):
        val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
        if val == 0:
            return 0  # collinear
        return 1 if val > 0 else 2  # 1 for clockwise, 2 for counterclockwise

    def left_hull(points):
        hull = []
        for p in points:
            while len(hull) >= 2 and orientation(hull[-2], hull[-1], p) != 2:
                hull.pop()
            hull.append(p)
        return hull

    def right_hull(points):
        hull = []
        for p in reversed(points):
            while len(hull) >= 2 and orientation(hull[-2], hull[-1], p) != 2:
                hull.pop()
            hull.append(p)
        return hull

    points.sort(key=lambda p: (p.x, p.y))
    upper_hull = left_hull(points)
    lower_hull = right_hull(points)

    return upper_hull[:-1] + lower_hull[:-1]

class ConvexHullApp:
    def __init__(self, master):
        self.master = master
        self.master.title("Convex Hull Algorithms")
        self.master.geometry("1000x800")

        self.algorithm_var = tk.StringVar()
        self.points = []

        algorithm_label = tk.Label(master, text="Select Algorithm:")
        algorithm_label.pack(pady=10)
        algorithms = ["Graham's Scan", "Quick Hull", "Jarvis March", "Brute Force", "Monotone Chain"]  
        algorithm_dropdown = ttk.Combobox(master, textvariable=self.algorithm_var, values=algorithms)
        algorithm_dropdown.pack(pady=10)
        algorithm_dropdown.set(algorithms[0])

        self.canvas = Canvas(master, bg="black", width=800, height=600)  
        self.canvas.pack(pady=20)

        calculate_button = tk.Button(master, text="Calculate Convex Hull", command=self.calculate_convex_hull)
        calculate_button.place(x=50, y=0, width=150, height=20)

        clear_button = tk.Button(master, text="Clear", command=self.clear_canvas, width=15, height=2) 
        clear_button.place(x=0, y=0, width=50, height=20)

        self.canvas.bind("<Button-1>", self.add_point)

    def add_point(self, event):
        x, y = event.x, event.y
        self.points.append(Point(x, y))
        self.canvas.create_oval(x - 3, y - 3, x + 3, y + 3, fill="white")

    def calculate_convex_hull(self):
        algorithm = self.algorithm_var.get()

        if not self.points:
            return

        if algorithm == "Graham's Scan":
            start_time = time.time()
            self.graham_scan()
            end_time = time.time()
            print(f"Graham's scan time complexity: O(n log n) - Elapsed Time: {end_time - start_time:.6f} seconds")
        elif algorithm == "Quick Hull":
            start_time = time.time()
            self.quick_hull()
            end_time = time.time()
            print(f"Quick Hull's time complexity: O(n log n) - Elapsed Time: {end_time - start_time:.6f} seconds")
        elif algorithm == "Jarvis March":
            start_time = time.time()
            self.jarvis_march()
            end_time = time.time()
            print(f"Jarvis March time complexity: O(nh) - Elapsed Time: {end_time - start_time:.6f} seconds")
        elif algorithm == "Brute Force":
            start_time = time.time()
            self.brute_force()
            end_time = time.time()
            print(f"Brute Force's time complexity: O(n log n) - Elapsed Time: {end_time - start_time:.6f} seconds")
        elif algorithm == "Monotone Chain":
            start_time = time.time()
            self.monotone_chain()
            end_time = time.time()
            print(f"Monotone Chain's time complexity: O(n log n) - Elapsed Time: {end_time - start_time:.6f} seconds")

    def clear_canvas(self):
        self.canvas.delete("all")
        self.points = []

    def graham_scan(self):
       
         
        global p0

        def next_to_top(S):
            return S[-2]

        def dist_sq(p1, p2):
            return ((p1.x - p2.x) * (p1.x - p2.x) +
                    (p1.y - p2.y) * (p1.y - p2.y))

        def orientation(p, q, r):
            val = ((q.y - p.y) * (r.x - q.x) -
                   (q.x - p.x) * (r.y - q.y))
            if val == 0:
                return 0  # collinear
            elif val > 0:
                return 1  # clock wise
            else:
                return 2  # counterclockwise

        def compare(p1, p2):
            o = orientation(p0, p1, p2)
            if o == 0:
                if dist_sq(p0, p2) >= dist_sq(p0, p1):
                    return -1
                else:
                    return 1
            else:
                if o == 2:
                    return -1
                else:
                    return 1

        n = len(self.points)

        ymin = self.points[0].y
        min_point = 0
        for i in range(1, n):
            y = self.points[i].y

            if (y < ymin) or (ymin == y and self.points[i].x < self.points[min_point].x):
                ymin = self.points[i].y
                min_point = i

        self.points[0], self.points[min_point] = self.points[min_point], self.points[0]

        p0 = self.points[0]
        self.points = sorted(self.points, key=cmp_to_key(compare))

        m = 1
        for i in range(1, n):
            while ((i < n - 1) and (orientation(p0, self.points[i], self.points[i + 1]) == 0)):
                i += 1

            self.points[m] = self.points[i]
            m += 1
        if m < 3:
            return

        S = []
        S.append(self.points[0])
        S.append(self.points[1])
        S.append(self.points[2])

        for i in range(3, m):
            while ((len(S) > 1) and
                   (orientation(next_to_top(S), S[-1], self.points[i]) != 2)):
                S.pop()
            S.append(self.points[i])

        convex_hull_points = []
        while S:
            p = S[-1]
            convex_hull_points.append((p.x, p.y))
            S.pop()

        for i in range(len(convex_hull_points) - 1):
            x1, y1 = convex_hull_points[i]
            x2, y2 = convex_hull_points[i + 1]
            self.canvas.create_line(x1, y1, x2, y2, fill="white", width=2)

        x1, y1 = convex_hull_points[-1]
        x2, y2 = convex_hull_points[0]
        self.canvas.create_line(x1, y1, x2, y2, fill="white", width=2)
       
        pass

    def quick_hull(self):
     
        global hull_quick

        def findSide(p1, p2, p):
            val = (p[1] - p1[1]) * (p2[0] - p1[0]) - (p2[1] - p1[1]) * (p[0] - p1[0])

            if val > 0:
                return 1
            if val < 0:
                return -1
            return 0

        def lineDist(p1, p2, p):
            return abs((p[1] - p1[1]) * (p2[0] - p1[0]) - (p2[1] - p1[1]) * (p[0] - p1[0]))

        def quickHullRec(a, n, p1, p2, side):
            ind = -1
            max_dist = 0

            for i in range(n):
                temp = lineDist(p1, p2, a[i])

                if (findSide(p1, p2, a[i]) == side) and (temp > max_dist):
                    ind = i
                    max_dist = temp

            if ind == -1:
                hull_quick.add((p1[0], p1[1]))
                hull_quick.add((p2[0], p2[1]))
                return

            quickHullRec(a, n, a[ind], p1, -findSide(a[ind], p1, p2))
            quickHullRec(a, n, a[ind], p2, -findSide(a[ind], p2, p1))

        point_tuples = [(point.x, point.y) for point in self.points]

        n = len(point_tuples)
        hull_quick = set()

        if n < 3:
            self.canvas.create_text(100, 100, text="Convex hull not possible", fill="white")
            return

        min_x = 0
        max_x = 0
        for i in range(1, n):
            if point_tuples[i][0] < point_tuples[min_x][0]:
                min_x = i
            if point_tuples[i][0] > point_tuples[max_x][0]:
                max_x = i

        quickHullRec(point_tuples, n, point_tuples[min_x], point_tuples[max_x], 1)
        quickHullRec(point_tuples, n, point_tuples[min_x], point_tuples[max_x], -1)

        for point in point_tuples:
            self.canvas.create_oval(point[0] - 3, point[1] - 3, point[0] + 3, point[1] + 3, fill="white")

        hull_points = list(hull_quick)
        hull_points.sort(key=lambda x: x[0])
        for i in range(len(hull_points) - 1):
            self.canvas.create_line(hull_points[i], hull_points[i + 1], fill="white")

        self.canvas.create_line(hull_points[-1], hull_points[0], fill="white")
       
        pass

    def jarvis_march(self):
        if len(self.points) < 3:
            return

        result_points = jarvis_march(self.points)
        result_coords = [(p.x, p.y) for p in result_points]

        for point in self.points:
            self.canvas.create_oval(point.x - 3, point.y - 3, point.x + 3, point.y + 3, fill="white")

        for i in range(len(result_coords) - 1):
            self.canvas.create_line(result_coords[i], result_coords[i + 1], fill="white")

        self.canvas.create_line(result_coords[-1], result_coords[0], fill="white")

        pass

    def brute_force(self):
      
        points = np.array([(point.x, point.y) for point in self.points])
        hull = ConvexHull(points)

        for simplex in hull.simplices:
            x1, y1 = points[simplex[0]]
            x2, y2 = points[simplex[1]]
            self.canvas.create_line(x1, y1, x2, y2, fill="white", width=2)
        pass

    def jarvis_march(self):
        convex_hull = jarvis_march(self.points)
        self.draw_convex_hull(convex_hull)

    def monotone_chain(self):
        convex_hull = monotone_chain(self.points)
        self.draw_convex_hull(convex_hull)

    def draw_convex_hull(self, convex_hull):
        self.canvas.delete("hull")
        for i in range(len(convex_hull)):
            start = convex_hull[i]
            end = convex_hull[(i + 1) % len(convex_hull)]
            self.canvas.create_line(start.x, start.y, end.x, end.y, fill="red", width=2, tags="hull")

if __name__ == "__main__":
    root = tk.Tk()
    app = ConvexHullApp(root)
    root.mainloop()