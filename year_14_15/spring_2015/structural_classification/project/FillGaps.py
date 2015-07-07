# -*- coding: UTF-8 -*-
__author__ = 'Tanya'


import numpy as np
import heapq

def distance (x1,x2):
    return sum(np.square(x1-x2))

def fill_by_mean(classes,data,mask):
    sz = data.shape
    n = sz[0] # observ
    m = sz[1] # attr
    for cl in classes:
        dt_cl = data[cl,:]
        msk_cl = [False for i in range(n)] # маска класса
        for i in cl:
            msk_cl[i] = True
        out_of_class = [x for x in range(n) if x not in cl]
        msk = mask  # маска пропусков в классе
        msk[out_of_class,:] = False
        for i in range(m):
            d_col=dt_cl[:,i]
            missed = np.where(msk[:,i])
            not_missed_in_class=np.where(np.logical_and(~ msk[:,i],msk_cl))
            data[missed,i]=np.mean(data[not_missed_in_class,i])
    return data

def get_new_partition(n,r,data):
    # import heapq
    # r - начальное количество кластеров
    # алгоритм "объединение" для заполнения пропусков
    r = r - 1
    ind = list(range(n))
    # print(ind, type(range(10)))
    x0 = data[0,:]
    ind.remove(0)
    new_ind = [0]
    new_ind_dist = [0]
    while len(ind)>1:
        min_dist = np.iinfo(np.int16).max
        min_dist_ind = -1
        # find min dist from x0
        for j in range(len(ind)):
            dist = distance(x0,data[ind[j],:])
            # print('dist', x0,data[ind[j],:])
            if dist<min_dist:
                min_dist = dist
                min_dist_ind = ind[j]
        # print(min_dist_ind)
        # if (min_dist_ind in ind):
            # print('ok ', len(ind),' ',min_dist,' ',dist)
        new_ind.append(min_dist_ind)
        new_ind_dist.append(min_dist)
        x0 = np.mean(data[new_ind],axis=0)
        ind.remove(min_dist_ind)
    new_ind.append(ind[0])
    new_ind_dist.append(distance(x0,data[ind[0],:]))
    # new_ind - упорядоченный "по сходсту" список индексов объектов
    # new_ind_dist - соответсвующие расстояния
    d_dist = np.diff(new_ind_dist)
    # точки, где максимальные "скачки" расстояния
    max_dist_points = np.sort(heapq.nlargest(r, range(len(d_dist)), d_dist.take))
    # далее разбиваем на классы, где границы содержатся в max_dist_points
    classes = [new_ind[0:max_dist_points[0]]]
    for i in range(len(max_dist_points)-1):
        classes.append(new_ind[max_dist_points[i]:max_dist_points[i+1]])
    classes.append(new_ind[max_dist_points[len(max_dist_points)-1]:len(new_ind)])
    #len(classes)
    # объединяем маленькие классы
    new_classes = []
    for i in range(len(classes)):
        if (len(classes[i])< 8):
            # смотрим, где меньше расстояние
            if i==0:
                new_classes.append(classes[0]+classes[1])
            elif i==len(classes)-1:
                 new_classes.append(classes[i-1]+classes[i])
            elif (d_dist[max_dist_points[i]]>d_dist[max_dist_points[i-1]]):
                new_classes.append(classes[i]+classes[i+1])
            else:
                 new_classes.append(classes[i]+classes[i-1])
        else:
            new_classes.append(classes[i])
    return new_classes

def fill_missing( data ):
    mask = np.isnan(data)
    sz = data.shape
    n = sz[0] # observ
    m = sz[1] # attr
    r = int(n/15)
    old_partition = []
    new_partition = [range(n)]
    while (old_partition!=new_partition):
        old_partition = new_partition
        # получаем начальное приближение пропусков без разеления на классы
        data = fill_by_mean(old_partition,data,mask)
        # выполняем класетризацию
        new_partition = get_new_partition(n,r,data)
    return data

