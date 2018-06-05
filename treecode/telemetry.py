# -*- coding: utf-8 -*-

import numpy as np
import treecode.energy_and_momentum as EM


def is_gravity_field_weak(particles, C_2):
    # Функция, выдающая ошибку, если гравитационное поле становится
    # слишком сильным для применения используемой модели
    ERROR_NAME = ''
    # Считаем величину phi / c^2, где phi - гравитационный потенциал
    array_phi = abs(particles[:, 10] / C_2)
    # если модуль phi / c^2 превышает определенное значение, то
    # гравитационное поле считаем сильным что выходит за границы
    # применимости используемой модели
    array_phi = array_phi >= 0.05
    if array_phi.any():
        ERROR_NAME = 'Strong gravity field error'
    return ERROR_NAME


def speed_limit(particles, C_2):
    # Функция, выдающая ошибку если скорость материальной
    # точки станет больше скорости света
    ERROR_NAME = ''
    v = np.zeros([np.size(particles, 0), 3])
    v = np.multiply(particles[:, 3:6], particles[:, 3:6])
    v_2 = v.sum(axis=1) >= C_2
    if v_2.any():
        ERROR_NAME = 'FTL error'
    return ERROR_NAME


def enegry_parameters(q, X):
    # Считаем количество вылетевших из системы частиц
    part_num = 0
    while X[part_num, 11] < 0:
        part_num += 1
        if part_num == (np.size(X, 0)):
            break
    momentum = EM.momentum_of_system(X)
    kinetic_energy = EM.system_kinetic_energy(X)
    potential_energy = EM.system_potential_energy(X)
    kinetic_in_volume = EM.kinetic_energy_in_volume(X, part_num)
    potential_in_volume = EM.potential_energy_in_volume(X, part_num)
    virial_coeff_system = - kinetic_energy / potential_energy
    virial_coeff_selected_volume = - kinetic_in_volume / potential_in_volume
#    momentum_in_volume = EM.momentum_of_system_in_volume(X, part_num)
    # Записываем энергию системы в отдельный массив
    # 1) номер шага
    # 2) кинетическая энергия всей системы
    # 3) потенциальная энергия всей системы
    # 4) полная энергия всей системы
    # 5) максимальная разница в кинетической энергии
    # между исследуемым шагом и предыдущим
    # 6) максимальная разница в потенциальной энергии
    # между исследуемым шагом и предыдущим
    # 7) импульс системы по оси X
    # 8) импульс системы по оси Y
    # 9) импульс системы по оси Z
    # 10) кинетическая энергия всех частиц в объеме
    # 11) потенциальная энергия всех частиц в объеме
    # 12) полная энергия всех частиц в объеме
    ENERGY = [q,
              kinetic_energy,
              potential_energy,
              EM.system_energy_Newton(X),
              EM.max_dT(X),
              EM.max_dU(X),
              momentum[0],
              momentum[1],
              momentum[2],
              kinetic_in_volume,
              potential_in_volume,
              EM.system_energy_in_volume(X, part_num),
              virial_coeff_system,
              virial_coeff_selected_volume]
    return ENERGY
