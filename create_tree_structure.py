# -*- coding: utf-8 -*-

import numpy as np
import math as m


def distribution(X0, n, d):
    # Распределение X_size частиц по ячейкам со стороной Distance
    # с последующей сортировкой по номерам ячеек (3.04.18)
    X_size = np.size(X0, 0)
    for N_local in range(X_size):
        n_x = int(m.floor(X0[N_local, 0] / d))
        n_y = int(m.floor(X0[N_local, 1] / d))
        n_z = int(m.floor(X0[N_local, 2] / d))
        if (n_x > n) or (n_y > n) or (n_z > n) or \
                (n_x < 0) or (n_y < 0) or (n_z < 0):
            X0[N_local, 11] = -1
        else:
            X0[N_local, 11] = n_x * n * n + n_y * n + n_z
    return X0[X0[:, 11].argsort(kind='mergesort')]


def particles_to_cell(Y, Y_size, order_n, n_max):
    # Функция, определяющая параметры самых малых ячеек из параметров
    # находящихся внутри частиц (13.04.18)
    R_local = np.zeros([n_max, 23])
    part_num = 0
    part_count = 0
    L_2 = 3 * Distance * Distance
    while Y[part_num, 11] < 0:
        part_num += 1
        if part_num == (np.size(Y, 0)):
            break
    for cell_num in range(n_max):
        R = np.zeros([12])
        if not part_num == Y_size:
            while Y[part_num, 11] == cell_num:
                R[0:3] += Y[part_num, 0:3] * Y[part_num, 6]
                R[3] += Y[part_num, 6]
                part_num += 1
                if part_num == Y_size:
                    break
        R[4] = part_count
        R[5] = part_num
        part_count = part_num
        d_xy = 0
        d_xz = 0
        d_yz = 0
        if not R[3] == 0:
            # Расчет положения центра масс ячейки
            R[0:3] = R[0:3] / R[3]
            # Расчет положения геометрического центра ячейки
            cell_x = cell_num // (n * n)
            R[6] = Distance * (0.5 + cell_x)
            R[7] = Distance * (0.5 + ((cell_num // n) - cell_x * n))
            R[8] = Distance * (0.5 + (cell_num % n))
            # Расчет квадрупольного момента для выбранной ячейки
            for s in range(int(R[4]), int(R[5])):
                R[9] += Y[s, 6] * (Y[s, 0] - R[0]) * (Y[s, 1] - R[1])
                R[10] += Y[s, 6] * (Y[s, 0] - R[0]) * (Y[s, 2] - R[2])
                R[11] += Y[s, 6] * (Y[s, 1] - R[1]) * (Y[s, 2] - R[2])
                d_xy += Y[s, 6] * Y[s, 0] * Y[s, 1]
                d_xz += Y[s, 6] * Y[s, 0] * Y[s, 2]
                d_yz += Y[s, 6] * Y[s, 1] * Y[s, 2]
            R[9:12] *= 3
            # Итоговый вид строки с параметрами ячейки
        R_local[cell_num] = [R[0], R[1], R[2], R[6], R[7], R[8],
                             R[3], R[9], R[10], R[11], L_2, order_n,
                             R[4], R[5], 0, 0, 0, 0, 0, 0,
                             d_xy, d_xz, d_yz]
    return R_local


def cells_to_cell(R_final, order_n, n_max):
    # Функция, вычисляющая параметры ячеек за счет
    # находящихся внутри ячеек с меньшим порядком (13.04.18)
    cell_length = Distance * (n / order_n)
    n_linear = order_n * 2
    n_total = int(m.pow(order_n, 3))
    R_local = np.zeros([n_total, 23])
    L_2 = 3 * Distance * Distance * n * n / (order_n * order_n)
    for cell_num in range(n_total):
        R = np.zeros([10])
        cell_x = cell_num // (order_n * order_n)
        cell_y = (cell_num // order_n) - cell_x * order_n
        cell_z = cell_num % order_n
        cell_num_0 = 2 * int(cell_x * n_linear * n_linear
                             + cell_y * n_linear + cell_z)
        Numbers = [cell_num_0, cell_num_0 + 1,
                   cell_num_0 + int(n_linear),
                   cell_num_0 + int(n_linear) + 1,
                   cell_num_0 + int(n_linear * n_linear),
                   cell_num_0 + int(n_linear * n_linear) + 1,
                   cell_num_0 + int(n_linear * n_linear + n_linear),
                   cell_num_0 + int(n_linear * n_linear + n_linear) + 1]
        d_xy = 0
        d_xz = 0
        d_yz = 0
#        D_xy = 0
#        D_xz = 0
#        D_yz = 0
        for u in range(8):
            # Определяем параметры центра масс
            R[0:3] += R_final[Numbers[u], 0:3] \
                    * R_final[Numbers[u], 6]
            R[3] += R_final[Numbers[u], 6]
            # Определяем доп. параметры, связанные с квадрупольным вкладом
#            D_xy += R_final[Numbers[u], 6]  \
#                * R_final[Numbers[u], 0] * R_final[Numbers[u], 1]
#            D_xz += R_final[Numbers[u], 6]  \
#                * R_final[Numbers[u], 0] * R_final[Numbers[u], 2]
#            D_yz += R_final[Numbers[u], 6]  \
#                * R_final[Numbers[u], 1] * R_final[Numbers[u], 2]
#            d_xy += R_final[Numbers[u], 20]
#            d_xz += R_final[Numbers[u], 21]
#            d_yz += R_final[Numbers[u], 22]
        if not R[3] == 0:
            # Расчет положения ЦМ и геометрического центра ячейки
            R[0:3] = R[0:3] / R[3]
            R[4] = cell_length * (0.5 + cell_x)
            R[5] = cell_length * (0.5 + cell_y)
            R[6] = cell_length * (0.5 + cell_z)
            # Расчет квадрупольного момента для выбранной ячейки
#            for s in range(8):
#                if not R_final[Numbers[s], 6] == 0:
#                    R[7] += R_final[Numbers[s], 6]          \
#                        * (R_final[Numbers[s], 0] - R[0])   \
#                        * (R_final[Numbers[s], 1] - R[1])
#                    R[8] += R_final[Numbers[s], 6]          \
#                        * (R_final[Numbers[s], 0] - R[0])   \
#                        * (R_final[Numbers[s], 2] - R[2])
#                    R[9] += R_final[Numbers[s], 6]          \
#                        * (R_final[Numbers[s], 1] - R[1])   \
#                        * (R_final[Numbers[s], 2] - R[2])
#            if (R[7] == 0) and (R[8] == 0) and (R[9] == 0):
#                R[7] = R_final[Numbers[:], 7].sum()
#                R[8] = R_final[Numbers[:], 8].sum()
#                R[9] = R_final[Numbers[:], 9].sum()
#            else:
#                R[7] += d_xy - D_xy
#                R[8] += d_xz - D_xz
#                R[9] += d_yz - D_yz
#                R[7:10] *= 3
#        Итоговый вид строки с параметрами ячейки
        R_local[cell_num] = [R[0], R[1], R[2], R[4], R[5], R[6], R[3],
                             R[7], R[8], R[9], L_2, order_n,
                             Numbers[0], Numbers[1], Numbers[2], Numbers[3],
                             Numbers[4], Numbers[5], Numbers[6], Numbers[7],
                             d_xy, d_xz, d_yz]
#    Корректируем номера "дочерних" ячеек
    R_local[:, 12:20] += n_total
    R_final[0:(-n_max), 12:20] += n_total
    return np.vstack((R_local, R_final))


def create_tree(Y, number_of_cells, cell_size):
    # Функция, позволяющая получить новые параметры частиц
    # из матрицы Y с помощью метода Tree code (13.04.18)
    global n
    global Distance
    Distance = cell_size
    order_n = number_of_cells
    n = number_of_cells
    Y_size = np.size(Y, 0)
    n_max = int(n * n * n)
    R_final = particles_to_cell(Y, Y_size, order_n, n_max)
    while order_n > 1:
        order_n *= 0.5
        R_final = cells_to_cell(R_final, order_n, n_max)
    return R_final
