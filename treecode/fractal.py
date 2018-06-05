# -*- coding: utf-8 -*-

import random
import numpy as np
import treecode.make_graph as graph


def sqr_dist(X, Y, n_x, n_y):
    dx = X[n_x, 0] - Y[n_y, 0]
    dy = X[n_x, 1] - Y[n_y, 1]
    dz = X[n_x, 2] - Y[n_y, 2]
    r2 = dx * dx + dy * dy + dz * dz
    return r2


def create_sphere(i, j, k,
                  INDENT,
                  CELL_LENGTH, NUMBER_OF_CELLS,
                  NUMBER_OF_SPHERES):
    LENGHT = CELL_LENGTH * NUMBER_OF_CELLS / NUMBER_OF_SPHERES
    SPHERE = np.zeros([5])
    x_r = LENGHT * (INDENT + i) / NUMBER_OF_SPHERES
    y_r = LENGHT * (INDENT + j) / NUMBER_OF_SPHERES
    z_r = LENGHT * (INDENT + k) / NUMBER_OF_SPHERES
    SPHERE[0] = x_r
    SPHERE[1] = y_r
    SPHERE[2] = z_r
    return SPHERE


def create_net_w_spheres(NUMBER_OF_CELLS, CELL_LENGTH, START_RADIUS):
    NUMBER_OF_SPHERES = NUMBER_OF_CELLS - 1
    INDENT = random.random() * CELL_LENGTH
    SPHERES = np.zeros([5, NUMBER_OF_SPHERES
                        * NUMBER_OF_SPHERES
                        * NUMBER_OF_SPHERES])
    NUMBER = 0
    for k in range(NUMBER_OF_SPHERES):
        for j in range(NUMBER_OF_SPHERES):
            for i in range(NUMBER_OF_SPHERES):
                SPHERES[NUMBER] = create_sphere(i, j, k,
                                                INDENT,
                                                CELL_LENGTH,
                                                NUMBER_OF_CELLS,
                                                NUMBER_OF_SPHERES)
                NUMBER += 1
    SPHERES[:, 3] = START_RADIUS
    return SPHERES


def make_one_step(X, SPHERES):
    X_SIZE = np.size(X, 0)
    NUMBER_OF_SPHERES = np.size(SPHERES, 0)
    SQ_RAD = SPHERES[0, 3] * SPHERES[0, 3]
    for n_s in range(NUMBER_OF_SPHERES):
        for n_x in range(X_SIZE):
            SQ_DISTANCE = sqr_dist(X, SPHERES, n_x, n_s)
            if SQ_DISTANCE <= SQ_RAD:
                SPHERES[n_s, 4] += X[n_x, 6]
    return SPHERES


def create_data(X, NUMBER_OF_CELLS, CELL_LENGTH, SMOOTH_DIST):
    START_RADIUS = 50 * SMOOTH_DIST
    NUMBER_OF_STEPS = int(1.5 * CELL_LENGTH / START_RADIUS) + 2
    SPHERES = create_net_w_spheres(NUMBER_OF_CELLS, CELL_LENGTH, START_RADIUS)
    NUMBER_OF_SPHERES = np.size(SPHERES, 0)
    M_PER_RADIUS = np.zeros([NUMBER_OF_STEPS, NUMBER_OF_SPHERES, 2])
    for n in range(NUMBER_OF_STEPS):
        SPHERES = make_one_step(X, SPHERES)
        M_PER_RADIUS[n] = SPHERES[:, 3:5]
        SPHERES[:, 3] += START_RADIUS
    return M_PER_RADIUS


def make_dependency(X, NUMBER_OF_CELLS, CELL_LENGTH, SMOOTH_DIST):
    M_PER_RADIUS = create_data(X, NUMBER_OF_CELLS, CELL_LENGTH, SMOOTH_DIST)
    graph.plot_dependency(M_PER_RADIUS)
