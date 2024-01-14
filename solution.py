from typing import List, Tuple

from data import runtests
import sys

import heapq as hp

sys.setrecursionlimit(100000)


class Node:
    def __init__(self, idx):
        self.idx = idx
        self.neighbours = {}
        self.is_articulation_point = False
        self.components = set()

    def connect_to(self, v, w):
        self.neighbours[v] = w


class Graph:
    def __init__(self, V: int):
        self.V = V
        self.adjacency_list: List[Node] = [None] + [Node(i) for i in range(1, V + 1)]

    def add_edge(self, a: int, b: int, w: int) -> None:
        self.adjacency_list[a].connect_to(b, w)


def find_articulation_points(graph: Graph):
    visited = [False] * (graph.V + 1)
    parent = [-1] * (graph.V + 1)
    low = [float("inf")] * (graph.V + 1)
    times = [float("inf")] * (graph.V + 1)
    articulation_points = set()
    biconnected_components_stack = []
    biconnected_components = []

    time = 0

    def DFSVisit(G: Graph, u: int):
        nonlocal time, parent, visited, articulation_points, low, biconnected_components_stack, times

        time += 1
        visited[u] = True
        times[u] = time
        low[u] = time

        children = 0
        for v in G.adjacency_list[u].neighbours:
            if not visited[v]:
                parent[v] = u
                biconnected_components_stack.append((u, v))
                DFSVisit(G, v)
                low[u] = min(low[u], low[v])
                children += 1
                if low[v] >= times[u] and parent[u] != -1:
                    articulation_points.add(u)
                    G.adjacency_list[u].is_articulation_point = True

                    biconnected_components.append([])
                    while biconnected_components_stack[-1] != (u, v):
                        biconnected_components[-1].append(biconnected_components_stack.pop())

                    biconnected_components[-1].append(biconnected_components_stack.pop())

                if children > 1 and parent[u] == -1:
                    biconnected_components.append([])
                    while biconnected_components_stack[-1] != (u, v):
                        biconnected_components[-1].append(biconnected_components_stack.pop())

                    biconnected_components[-1].append(biconnected_components_stack.pop())

            elif parent[u] != v and times[v] < times[u]:
                biconnected_components_stack.append((u, v))
                low[u] = times[v]
            elif parent[u] != v:
                low[u] = min(low[u], times[v])

        if parent[u] == -1 and children > 1:
            articulation_points.add(u)
            G.adjacency_list[u].is_articulation_point = True

    for u in range(1, graph.V + 1):
        if not visited[u]:
            DFSVisit(graph, u)

            biconnected_components.append([])
            while biconnected_components_stack:
                biconnected_components[-1].append(biconnected_components_stack.pop())

    return articulation_points, biconnected_components


class Metric:
    def __init__(self, distance, arches):
        self.distance = distance
        self.arches = arches

    def __lt__(self, other):
        if self.arches > other.arches:
            return True
        elif self.arches < other.arches:
            return False

        return self.distance < other.distance


def dijkstra(G: Graph, s: int):
    metrics = [Metric(float("inf"), float("-inf")) for _ in range(G.V + 1)]
    metrics[s] = Metric(0, 0)
    parents: List[int] = [-1000000000]*(G.V + 1)
    visited = [False] * (G.V + 1)

    q: List[(Metric, int)] = []

    hp.heappush(q, (Metric(0, 0), s))

    while len(q) != 0:
        metric, v = hp.heappop(q)
        visited[v] = True

        for u, weight in G.adjacency_list[v].neighbours.items():
            if visited[u]:
                continue

            if v != s and G.adjacency_list[v].is_articulation_point:
                if not G.adjacency_list[u].components & G.adjacency_list[parents[v]].components:
                    if Metric(metric.distance + weight, metric.arches + 1) < metrics[u]:
                        metrics[u] = Metric(metric.distance + weight, metric.arches + 1)
                        parents[u] = v
                        hp.heappush(q, (metrics[u], u))
                        continue
            if Metric(metric.distance + weight, metric.arches) < metrics[u]:
                metrics[u] = Metric(metric.distance + weight, metric.arches)
                parents[u] = v
                hp.heappush(q, (metrics[u], u))

    return metrics, parents


def parade(N: int, streets: List[Tuple[int, int, int]]):
    """
    N: liczba skrzyżowań - V lista wierzchołków w grafie
    streets: lista ulic, każda ulica to krotka (a, b, t) - lista krawędzi w grafie
    """
    graph = Graph(N)
    for a, b, t in streets:
        graph.add_edge(a, b, t)
        graph.add_edge(b, a, t)

    articulation_points, biconnected_components = find_articulation_points(graph)
    for i, component in enumerate(biconnected_components):
        for u, v in component:
            graph.adjacency_list[u].components.add(i)
            graph.adjacency_list[v].components.add(i)

    global_best = Metric(float("inf"), float("-inf"))

    for articulation_point in articulation_points:

        for neighbour in graph.adjacency_list[articulation_point].neighbours:
            distances, parents = dijkstra(graph, neighbour)

            for i, metric in enumerate(distances):
                if metric < global_best:
                    global_best = metric

    return global_best.arches, global_best.distance


runtests(parade)
