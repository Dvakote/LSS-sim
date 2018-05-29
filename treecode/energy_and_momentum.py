# -*- coding: utf-8 -*-

import numpy as np


def momentum_of_system(Y):
    # Функция, определяющая импульс всей системы и выводящая его в строку
    P = np.zeros([np.size(Y, 0), 3])
    P[:, 0] = np.multiply(Y[:, 3], Y[:, 6])
    P[:, 1] = np.multiply(Y[:, 4], Y[:, 6])
    P[:, 2] = np.multiply(Y[:, 5], Y[:, 6])
    momentum = P.sum(axis=0)
    return momentum


def momentum_of_particles(Y):
    # Функция, определяющая импульс всех материальных точек
    P = np.zeros([np.size(Y, 0), 3])
    P[:, 0] = np.multiply(Y[:, 3], Y[:, 6])
    P[:, 1] = np.multiply(Y[:, 4], Y[:, 6])
    P[:, 2] = np.multiply(Y[:, 5], Y[:, 6])
    if np.size(Y, 0) > 10:
        print('Импульсы всех материальных точек сохранены в файл')
        np.savetxt('Импульсы материальных точек.txt', P)
    else:
        print(P)


def set_system_momentum_to_0(Y):
    size_Y = np.size(Y, 0)
    speed_change = np.zeros([size_Y, 3])
    total_momentum = momentum_of_system(Y)
    total_momentum = total_momentum / size_Y
    for i in range(size_Y):
        speed_change[i] = total_momentum
    speed_change[:, 0] /= Y[:, 6]
    speed_change[:, 1] /= Y[:, 6]
    speed_change[:, 2] /= Y[:, 6]
    Y[:, 3:6] = Y[:, 3:6] - speed_change
    return Y


def kinetic_energy_Newton(Y):
    # Функция, определяющая кинетическую энергию каждой частицы
    V = np.multiply(Y[:, 3:6], Y[:, 3:6])
    E = V.sum(axis=1)
    E = np.multiply(E[:], Y[:, 6])
    E *= 0.5
    return E


def potential_energy_Newton(Y):
    # Функция, определяющая потенциальную энергию каждой частицы
    E = np.multiply(Y[:, 10], Y[:, 6])
    E *= 0.5
    return E


def system_kinetic_energy(Y):
    # Функция, определяющая полную энергию системы
    E = kinetic_energy_Newton(Y)
    E = E.sum(axis=0)
    return E


def system_potential_energy(Y):
    E = potential_energy_Newton(Y)
    E = E.sum(axis=0)
    return E


def system_energy_Newton(Y):
    # Функция, определяющая полную энергию системы
    E = system_kinetic_energy(Y)
    E += system_potential_energy(Y)
    return E


def max_dT(Y):
    # Функция, определяющая максимальную разницу
    # кинетической энергии частиц за шаг
    E = kinetic_energy_Newton(Y)
    E = E - Y[:, 12]
    dE_plus = np.amax(E)
    dE_minus = np.amin(E)
    if abs(dE_minus) > dE_plus:
        dE = dE_minus
    else:
        dE = dE_plus
    return dE


def max_dU(Y):
    # Функция, определяющая максимальную разницу
    # потенциальной энергии частиц за шаг
    E = potential_energy_Newton(Y)
    E = E - Y[:, 13]
    dE_plus = np.amax(E)
    dE_minus = np.amin(E)
    if abs(dE_minus) > dE_plus:
        dE = dE_minus
    else:
        dE = dE_plus
    return dE
