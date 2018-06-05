# -*- coding: utf-8 -*-
"""
Редактор Spyder

@author: Дмитрий Мелкозеров
"""

import numpy as np
import tree_code_math as tcm


def N_body_direct(X0, SYSTEM_PARAMETERS):
    # Ньютоновская гравитация, метод частица-частица
    smooth = float(SYSTEM_PARAMETERS[1])
    GRAVITY_TYPE = str(SYSTEM_PARAMETERS[2])
    ALPHA = float(SYSTEM_PARAMETERS[3])
    LAMBDA = float(SYSTEM_PARAMETERS[4])
    GAMMA = float(SYSTEM_PARAMETERS[5])
    total_part = np.size(X0, 0)
    A = np.zeros((total_part, 4))
    part_num = 0
    while X0[part_num, 11] < 0:
        part_num += 1
        if part_num == total_part:
            break
    if GRAVITY_TYPE == 'Newton':
        A = tcm.g_force_Newton(X0, part_num, total_part, smooth)
    elif GRAVITY_TYPE == 'plusYukawa':
        A = tcm.g_force_plusYukawa(X0, part_num, total_part,
                                   smooth, ALPHA, LAMBDA)
    elif GRAVITY_TYPE == 'BD':
        A = tcm.g_force_BD(X0, part_num, total_part, smooth, GAMMA)
    # for l in range(total_part):
    #   if not part_num == l:
    #       A[l] = g_force_Newton(X0, part_num, l, smooth)
    X0[:, 7:11] += A
    return X0


def branch_to_leafes(Mass_center, current_cell, cell_num, Numbers):
    # Функция, рассчитывающая гравитационное воздействие на частицы в
    # ячейке current_cell со стороны ячейки cell_num
    if Mass_center[current_cell, 9] == n:
        Numbers.append(current_cell)
    else:
        if not Mass_center[int(Mass_center[current_cell, 10]), 6] == 0:
            Numbers = branch_to_leafes(Mass_center,
                                       int(Mass_center[current_cell, 10]),
                                       cell_num, Numbers)
        if not Mass_center[int(Mass_center[current_cell, 11]), 6] == 0:
            Numbers = branch_to_leafes(Mass_center,
                                       int(Mass_center[current_cell, 11]),
                                       cell_num, Numbers)
        if not Mass_center[int(Mass_center[current_cell, 12]), 6] == 0:
            Numbers = branch_to_leafes(Mass_center,
                                       int(Mass_center[current_cell, 12]),
                                       cell_num, Numbers)
        if not Mass_center[int(Mass_center[current_cell, 13]), 6] == 0:
            Numbers = branch_to_leafes(Mass_center,
                                       int(Mass_center[current_cell, 13]),
                                       cell_num, Numbers)
        if not Mass_center[int(Mass_center[current_cell, 14]), 6] == 0:
            Numbers = branch_to_leafes(Mass_center,
                                       int(Mass_center[current_cell, 14]),
                                       cell_num, Numbers)
        if not Mass_center[int(Mass_center[current_cell, 15]), 6] == 0:
            Numbers = branch_to_leafes(Mass_center,
                                       int(Mass_center[current_cell, 15]),
                                       cell_num, Numbers)
        if not Mass_center[int(Mass_center[current_cell, 16]), 6] == 0:
            Numbers = branch_to_leafes(Mass_center,
                                       int(Mass_center[current_cell, 16]),
                                       cell_num, Numbers)
        if not Mass_center[int(Mass_center[current_cell, 17]), 6] == 0:
            Numbers = branch_to_leafes(Mass_center,
                                       int(Mass_center[current_cell, 17]),
                                       cell_num, Numbers)
    return Numbers


def tree_branch(Particles, Mass_center, current_cell, cell_num, A):
    # Функция, определяющая дальнейший алгоритм действий, исходя из
    # заданного критерия раскрытия ячеек (15.04.18)
    sqr_dist = tcm.sqr_dist(Mass_center, current_cell, cell_num)
    if sqr_dist > Mass_center[cell_num, 8]:
        Numbers_of_cells = []
        Numbers_of_cells = branch_to_leafes(Mass_center,
                                            current_cell, cell_num,
                                            Numbers_of_cells)
        Numbers_of_particles = []
        for l in Numbers_of_cells:
            for k in range(int(Mass_center[l, 10]),
                           int(Mass_center[l, 11])):
                Numbers_of_particles.append(k)
        if GRAVITY_TYPE == 'Newton':
            for p in Numbers_of_particles:
                A[p] += tcm.int_C_to_P(Particles, Mass_center,
                                       p, cell_num)
        elif GRAVITY_TYPE == 'plusYukawa':
            print('Yukawa')
            for p in Numbers_of_particles:
                A[p] += tcm.interaction_C2P_plus_Yukawa(Particles,
                                                        Mass_center,
                                                        alpha, lambd,
                                                        p, cell_num)
        elif GRAVITY_TYPE == 'BD':
            for p in Numbers_of_particles:
                A[p] += tcm.interaction_C2P_Brance_Dicke(Particles,
                                                         Mass_center,
                                                         gamma,
                                                         p, cell_num)
#            A[p] += int_C_to_P(Particles, Mass_center, p, cell_num)
    else:
        if Mass_center[cell_num, 9] == n:
            n1 = int(Mass_center[current_cell, 10])
            n2 = int(Mass_center[current_cell, 11])
            n_1 = int(Mass_center[cell_num, 10])
            n_2 = int(Mass_center[cell_num, 11])
            if GRAVITY_TYPE == 'Newton':
                A[n1:n2] += tcm.int_Ps_to_P(Particles, n1, n2,
                                            n_1, n_2, eps_smooth)
            elif GRAVITY_TYPE == 'plusYukawa':
                A[n1:n2] += tcm.interaction_Ps2P_plus_Yukawa(Particles,
                                                             n1, n2,
                                                             n_1, n_2,
                                                             eps_smooth,
                                                             alpha, lambd)
            elif GRAVITY_TYPE == 'BD':
                A[n1:n2] += tcm.interaction_Ps2P_Brance_Dicke(Particles,
                                                              n1, n2,
                                                              n_1, n_2,
                                                              eps_smooth,
                                                              gamma)
#            for Part_num in range(n1, n2):
#                A[Part_num] += int_Ps_to_P(Particles, Part_num, n_1, n_2)
        else:
            A = main_tree_branch(Particles, Mass_center,
                                 current_cell, cell_num, A)
    return A


def sub_tree_branch(Particles, Mass_center, current_cell, cell_num, A):
    # Функция, которая задает ячейки, воздействующие на частицы (15.04.18)
    if not Mass_center[int(Mass_center[cell_num, 10]), 6] == 0:
        A = tree_branch(Particles, Mass_center, current_cell,
                        int(Mass_center[cell_num, 10]), A)
    if not Mass_center[int(Mass_center[cell_num, 11]), 6] == 0:
        A = tree_branch(Particles, Mass_center, current_cell,
                        int(Mass_center[cell_num, 11]), A)
    if not Mass_center[int(Mass_center[cell_num, 12]), 6] == 0:
        A = tree_branch(Particles, Mass_center, current_cell,
                        int(Mass_center[cell_num, 12]), A)
    if not Mass_center[int(Mass_center[cell_num, 13]), 6] == 0:
        A = tree_branch(Particles, Mass_center, current_cell,
                        int(Mass_center[cell_num, 13]), A)
    if not Mass_center[int(Mass_center[cell_num, 14]), 6] == 0:
        A = tree_branch(Particles, Mass_center, current_cell,
                        int(Mass_center[cell_num, 14]), A)
    if not Mass_center[int(Mass_center[cell_num, 15]), 6] == 0:
        A = tree_branch(Particles, Mass_center, current_cell,
                        int(Mass_center[cell_num, 15]), A)
    if not Mass_center[int(Mass_center[cell_num, 16]), 6] == 0:
        A = tree_branch(Particles, Mass_center, current_cell,
                        int(Mass_center[cell_num, 16]), A)
    if not Mass_center[int(Mass_center[cell_num, 17]), 6] == 0:
        A = tree_branch(Particles, Mass_center, current_cell,
                        int(Mass_center[cell_num, 17]), A)
    return A


def main_tree_branch(Particles, Mass_center, current_cell, cell_num, A):
    # Функция, задающая ячейки, частицы в которых
    # будут испытывать воздействие (15.04.18)
    if not Mass_center[int(Mass_center[current_cell, 10]), 6] == 0:
        A = sub_tree_branch(Particles, Mass_center,
                            int(Mass_center[current_cell, 10]), cell_num, A)
    if not Mass_center[int(Mass_center[current_cell, 11]), 6] == 0:
        A = sub_tree_branch(Particles, Mass_center,
                            int(Mass_center[current_cell, 11]), cell_num, A)
    if not Mass_center[int(Mass_center[current_cell, 12]), 6] == 0:
        A = sub_tree_branch(Particles, Mass_center,
                            int(Mass_center[current_cell, 12]), cell_num, A)
    if not Mass_center[int(Mass_center[current_cell, 13]), 6] == 0:
        A = sub_tree_branch(Particles, Mass_center,
                            int(Mass_center[current_cell, 13]), cell_num, A)
    if not Mass_center[int(Mass_center[current_cell, 14]), 6] == 0:
        A = sub_tree_branch(Particles, Mass_center,
                            int(Mass_center[current_cell, 14]), cell_num, A)
    if not Mass_center[int(Mass_center[current_cell, 15]), 6] == 0:
        A = sub_tree_branch(Particles, Mass_center,
                            int(Mass_center[current_cell, 15]), cell_num, A)
    if not Mass_center[int(Mass_center[current_cell, 16]), 6] == 0:
        A = sub_tree_branch(Particles, Mass_center,
                            int(Mass_center[current_cell, 16]), cell_num, A)
    if not Mass_center[int(Mass_center[current_cell, 17]), 6] == 0:
        A = sub_tree_branch(Particles, Mass_center,
                            int(Mass_center[current_cell, 17]), cell_num, A)
    return A


def begin_tree(Particles, Mass_center,
               current_cell, SYSTEM_PARAMETERS):
    global n
    global eps_smooth
    global GRAVITY_TYPE
    global alpha
    global lambd
    global gamma
    n = int(SYSTEM_PARAMETERS[0])
    eps_smooth = float(SYSTEM_PARAMETERS[1])
    GRAVITY_TYPE = str(SYSTEM_PARAMETERS[2])
    alpha = float(SYSTEM_PARAMETERS[3])
    lambd = float(SYSTEM_PARAMETERS[4])
    gamma = float(SYSTEM_PARAMETERS[5])
    A = np.zeros([np.size(Particles, 0), 4])
    A = sub_tree_branch(Particles, Mass_center, current_cell, 0, A)
    return A


# def part_distance(Particle_1, Particle_2, Number_1, Nubmer_2):
#    # Функция, которая выдает растояние между частицам 1 и 2
#    delta_x = Particle_1[Number_1, 0] - Particle_2[Nubmer_2, 0]
#    delta_y = Particle_1[Number_1, 1] - Particle_2[Nubmer_2, 1]
#    delta_z = Particle_1[Number_1, 2] - Particle_2[Nubmer_2, 2]
#    return m.sqrt(delta_x * delta_x
#                  + delta_y * delta_y
#                  + delta_z * delta_z)
#
#
# def smooth_distance(Particles, Number_1, Nubmer_2):
#    # Функция, выдающая растояние между частицам 1 и 2
#    delta_x = Particles[Number_1, 0] - Particles[Nubmer_2, 0]
#    delta_y = Particles[Number_1, 1] - Particles[Nubmer_2, 1]
#    delta_z = Particles[Number_1, 2] - Particles[Nubmer_2, 2]
#    delta_soft = m.sqrt(delta_x * delta_x + delta_y * delta_y
#                        + delta_z * delta_z + eps_smooth * eps_smooth)
#    return delta_soft  # * delta_soft * delta_soft
#
#
# def g_force_Newton(Particles, part_num, l, smooth):
#    # Ускорение по Ньютону
#    a = np.zeros([4])
#    r_1 = tcm.smooth_distance(Particles, part_num, l, smooth)
# #    smooth_distance(Particles, part_num, l)
#    r_3 = r_1 * r_1 * r_1
#    a[0] = Particles[part_num, 6] \
#        * (Particles[part_num, 0] - Particles[l, 0]) / r_3
#    a[1] = Particles[part_num, 6] \
#        * (Particles[part_num, 1] - Particles[l, 1]) / r_3
#    a[2] = Particles[part_num, 6] \
#        * (Particles[part_num, 2] - Particles[l, 2]) / r_3
#    a[3] = - Particles[part_num, 6] / r_1
#    return a
#
#
# def quadrupole(Mass_center, num, r_1, r_3, delta_x, delta_y, delta_z):
#    # Функция, расчитывающая квадрупольный вклад
#    r_5 = r_3 * r_1 * r_1
#    r_7 = r_5 * r_1 * r_1
#    DR = (Mass_center[num, 7] * delta_x * delta_y
#          + Mass_center[num, 8] * delta_x * delta_z
#          + Mass_center[num, 9] * delta_y * delta_z) * 5
#    a_x = - (Mass_center[num, 7] * delta_y + Mass_center[num, 8] * delta_z) \
#        / r_5 + DR * delta_x / r_7
#    a_y = - (Mass_center[num, 7] * delta_x + Mass_center[num, 9] * delta_z) \
#        / r_5 + DR * delta_y / r_7
#    a_z = - (Mass_center[num, 8] * delta_x + Mass_center[num, 9] * delta_y) \
#        / r_5 + DR * delta_z / r_7
#    phi = DR / (5 * r_5)
#    return np.array([a_x, a_y, a_z, - phi])
#
#
# def int_C_to_P(Particles, Mass_center, Part_num, cell_num):
#    # Функция, рассчитывающая ускорение частицы под номером Part_num,
#    # полученное за счет гравитационного мультипольного взаимодействия с
#    # частицами в ячейке с номером cell_num.
#    r_1 = part_distance(Particles, Mass_center, Part_num, cell_num)
#    r_3 = r_1 * r_1 * r_1
#    delta_x = Mass_center[cell_num, 0] - Particles[Part_num, 0]
#    delta_y = Mass_center[cell_num, 1] - Particles[Part_num, 1]
#    delta_z = Mass_center[cell_num, 2] - Particles[Part_num, 2]
#    cell_to_body = np.array([delta_x, delta_y, delta_z, 0])
#    cell_to_body[0:3] *= Mass_center[cell_num, 6] / r_3
#    cell_to_body[3] = - Mass_center[cell_num, 6] / r_1
# #    cell_to_body += quadrupole(Mass_center, cell_num, r_1, r_3,
# #                               delta_x, delta_y, delta_z)
#    return cell_to_body
#
#
# def int_Ps_to_P(Particles, Part_num, n1, n2):
#    # Функция, рассчитывающая ускорение частицы под номером Part_num,
#    # полученное за счет гравитационного взаимодействия с частицами
#    # в ячейке с номером cell_num. (Для использования в методе  Tree code)
#    a_x = 0
#    a_y = 0
#    a_z = 0
#    phi = 0
#    if (Part_num >= n1) and (Part_num < n2):
#        for num in range(n1, n2):
#            if not num == Part_num:
#                r_1 = smooth_distance(Particles, Part_num, num)
#                r_3 = r_1 * r_1 * r_1
#                a_x += Particles[num, 6] \
#                    * (Particles[num, 0] - Particles[Part_num, 0]) / r_3
#                a_y += Particles[num, 6] \
#                    * (Particles[num, 1] - Particles[Part_num, 1]) / r_3
#                a_z += Particles[num, 6] \
#                    * (Particles[num, 2] - Particles[Part_num, 2]) / r_3
#                phi += Particles[num, 6] / r_1
#    else:
#        for num in range(n1, n2):
#            r_1 = smooth_distance(Particles, Part_num, num)
#            r_3 = r_1 * r_1 * r_1
#            a_x += Particles[num, 6] \
#                * (Particles[num, 0] - Particles[Part_num, 0]) / r_3
#            a_y += Particles[num, 6] \
#                * (Particles[num, 1] - Particles[Part_num, 1]) / r_3
#            a_z += Particles[num, 6] \
#                * (Particles[num, 2] - Particles[Part_num, 2]) / r_3
#            phi += Particles[num, 6] / r_1
#    return np.array([a_x, a_y, a_z, - phi])
