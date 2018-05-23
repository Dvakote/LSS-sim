# -*- coding: utf-8 -*-

import numpy as np


def momentum_of_system(Y):
    # Функция, определяющая импульс всей системы и выводящая его в строку
    P = np.zeros([np.size(Y, 0), 3])
    P[:, 0] = np.multiply(Y[:, 3], Y[:, 6])
    P[:, 1] = np.multiply(Y[:, 4], Y[:, 6])
    P[:, 2] = np.multiply(Y[:, 5], Y[:, 6])
    print('Полный импульс системы ', P.sum(axis=0))


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


def kinetic_energy_Newton(Y):
    # Функция, определяющая кинетическую энергию каждой частицы
    V = np.multiply(Y[:, 3:6], Y[:, 3:6])
    E = V.sum(axis=1)
    E = np.multiply(E[:], Y[:, 6])
    E /= 2
    return E


def potential_energy_Newton(Y):
    # Функция, определяющая кинетическую энергию каждой частицы
    E = np.multiply(Y[:, 10], Y[:, 6])
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
    E = E + system_potential_energy(Y)
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
