# -*- coding: utf-8 -*-

import random as r
import numpy as np


def parameters_test(h, p, l,
                    indent_i, indent_j,
                    indent_k, period,
                    i_test, j_test, k_test):
    # Подфункция, позволяющая сгенерировать определенные
    # параметры для тела
    x = Distance * (indent_i + h * period) / i_test
    y = Distance * (indent_j + p * period) / j_test
    z = Distance * (indent_k + l * period) / k_test
    # Распределение скоростей и масс считаем нормальным
    Vx = r.normalvariate(0, 4) * v_avg
    Vy = r.normalvariate(0, 4) * v_avg
    Vz = r.normalvariate(0, 4) * v_avg
    mass = abs(m_avg)
    Sum = np.array([x, y, z, Vx, Vy, Vz, mass, 0, 0, 0, 0, 0, 0, 0])
    return Sum


def randomize_parameters():
    # Подфункция, позволяющая сгенерировать случайные параметры для тела
    # внутри рассматриваемого объема
    x = r.random() * n * Distance
    y = r.random() * n * Distance
    z = r.random() * n * Distance
#   Распределение скоростей и масс считаем нормальным
#   (пока что квадратичное отклонение выбрано наугад)
    Vx = r.normalvariate(0, 4) * v_avg
    Vy = r.normalvariate(0, 4) * v_avg
    Vz = r.normalvariate(0, 4) * v_avg
    mass = abs(r.normalvariate(m_avg, 0.5*m_avg))
    Sum = np.array([x, y, z, Vx, Vy, Vz, mass, 0, 0, 0, 0, 0, 0, 0])
    return Sum


def randomize_ellipsoid(a_inp, b_inp, c_inp, w_x, w_y, w_z):
    # Подфункция, позволяющая сгенерировать случайные параметры для тела
    # внутри эллипсоида заданного объема
    x_r = 0
    y_r = 0
    z_r = 0
    particle_not_generated = True
    while particle_not_generated:
        x_r = r.random()
        y_r = r.random()
        z_r = r.random()
        x_el = (2 * x_r - 1) / a_inp
        y_el = (2 * y_r - 1) / b_inp
        z_el = (2 * z_r - 1) / c_inp
        ellipsoid = x_el * x_el + y_el * y_el + z_el * z_el
        if ellipsoid <= 1:
            particle_not_generated = False
    center = n * Distance / 2
    x = (x_r + 0.5) * center
    y = (y_r + 0.5) * center
    z = (z_r + 0.5) * center
    d_x = x - center
    d_y = y - center
    d_z = z - center
#   Распределение скоростей и масс считаем нормальным
#   (пока что квадратичное отклонение выбрано наугад)
    Vx = r.normalvariate(0, 3) * v_avg + w_y * d_z - w_z * d_y
    Vy = r.normalvariate(0, 3) * v_avg + w_z * d_x - w_x * d_z
    Vz = r.normalvariate(0, 3) * v_avg + w_x * d_y - w_y * d_x
    mass = abs(r.normalvariate(m_avg, 0.5*m_avg))
    Sum = np.array([x, y, z, Vx, Vy, Vz, mass, 0, 0, 0, 0, 0, 0, 0])
    return Sum


def birth_test(inp_parametrs):
    # Функция, создающая i*j*k тел
    # Сначала создаем массив нулей, а затем заполняем его;
    # тела находятся по первому индексу, параметры - по второму
    i_test = int(inp_parametrs[11])
    j_test = int(inp_parametrs[12])
    k_test = int(inp_parametrs[13])
    indent_i = inp_parametrs[14]
    indent_j = inp_parametrs[15]
    indent_k = inp_parametrs[16]
    period = inp_parametrs[17]
    test_particles = np.zeros((i_test * j_test * k_test, 14))
    Num = 0
    for l in range(k_test):
        for p in range(j_test):
            for h in range(i_test):
                test_particles[Num] = parameters_test(h, p, l,
                                                      indent_i, indent_j,
                                                      indent_k, period,
                                                      i_test, j_test, k_test)
                Num += 1
    return test_particles


def birth_random(inp_parametrs):
    # Функция, создающая "body_count" тел
    # Сначала создаем массив нулей, а затем заполняем его;
    # тела находятся по первому индексу, параметры - по второму
    random_particles = np.zeros((int(inp_parametrs[0]), 14))
    for l in range(inp_parametrs[0]):
        random_particles[l] = randomize_parameters()
    return random_particles


def birth_ellipsoid(inp_parametrs):
    # Функция, создающая "body_count" тел
    # Сначала создаем массив нулей, а затем заполняем его;
    # тела находятся по первому индексу, параметры - по второму
    a_inp = inp_parametrs[5]
    b_inp = inp_parametrs[6]
    c_inp = inp_parametrs[7]
    w_x = inp_parametrs[8]
    w_y = inp_parametrs[9]
    w_z = inp_parametrs[10]
    random_particles = np.zeros([int(inp_parametrs[0]), 14])
    for l in range(inp_parametrs[0]):
        random_particles[l] = randomize_ellipsoid(a_inp, b_inp, c_inp,
                                                  w_x, w_y, w_z)
    return random_particles


def generate_system(config_name, inp_parametrs):
    # Вспомогательная функция, необходимая для работы с основной
    # частью программы
    global m_avg
    global v_avg
    global n
    global Distance
    m_avg = inp_parametrs[1]
    v_avg = inp_parametrs[2]
    n = int(inp_parametrs[3])
    Distance = inp_parametrs[4]
    if config_name == 'cube':
        empty_config = birth_test(inp_parametrs)
    elif config_name == 'random':
        empty_config = birth_random(inp_parametrs)
    elif config_name == 'ellipsoid':
        if (inp_parametrs[5] == 0) or (inp_parametrs[6] == 0) \
                or (inp_parametrs[7] == 0):
            empty_config = 'zero volume'
        else:
            empty_config = birth_ellipsoid(inp_parametrs)
    return empty_config
