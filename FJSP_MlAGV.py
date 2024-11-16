"""
需要完成的编程工作：
1、初始化规则，5分
2、VNS规则，0分
3、节能规则，0分
4、正交实验，0分这部分也是前面弄好差不多就好了
5、MILP模型验证：目前少能耗部分的编码，其实大概也编写完成了，7分
6、跟其他五种算法的对比，需要编写其他五种算法的代码，顺便可以整理代码的框架，0分
7、AGV的载重以及AGV的数量的实验对比：如果前面的代码编好，这部分就很容易了，0分
以上的工作其实就是两部分工作，一部分把自己本身的代码写完成，另一部分就是编写其他算法的代码
"""
import random
import copy
from Data_transfer import Get_Problem
from FJSP_agv import FJSP_agv,Job,Machine,Agv
from itertools import chain
import numpy as np
import math
import sys
import traceback
from myplot import save_plot_main_data
import Logging
logger = Logging.getLogger(__name__)

# 快速非支配排序
def fast_non_dominated_sort(values1, values2):
    S = [[] for i in range(0, len(values1))]
    front = [[]]
    n = [0 for i in range(0, len(values1))]
    rank = [0 for i in range(0, len(values1))]

    for p in range(0, len(values1)):
        S[p] = []
        n[p] = 0
        for q in range(0, len(values1)):
            if (values1[p] <= values1[q] and values2[p] <= values2[q]):
                if q not in S[p]:
                    S[p].append(q)
            elif (values1[q] <= values1[p] and values2[q] <= values2[p]):
                n[p] = n[p] + 1
        if n[p] == 0:
            rank[p] = 0
            if p not in front[0]:
                front[0].append(p)

    i = 0
    while (front[i] != []):
        Q = []
        for p in front[i]:
            for q in S[p]:
                n[q] = n[q] - 1
                if (n[q] == 0):
                    rank[q] = i + 1
                    if q not in Q:
                        Q.append(q)
        i = i + 1
        front.append(Q)

    del front[len(front) - 1]
    return front


# 拥挤度选择
def crowding_distance(values1, values2, I):
    distance = [0 for i in range(0, len(I))]
    # sorted1 = sort_by_values(I, values1[:])
    # sorted2 = sort_by_values(I, values2[:])
    sorted1 = sorted(I, key=lambda x: values1[x])  
    sorted2 = sorted(I, key=lambda x: values2[x]) 
    distance[0] = 4444444444444444
    distance[len(I) - 1] = 4444444444444444
    for k in range(1, len(I) - 1):
        distance[k] = distance[k] + (values1[sorted1[k + 1]] - values1[sorted1[k - 1]]) / (
                max(values1) - min(values1))
    for k in range(1, len(I) - 1):
        distance[k] = distance[k] + (values2[sorted2[k + 1]] - values2[sorted2[k - 1]]) / (
                max(values2) - min(values2))#这边我进行了修整，不知道对不对，跟参考的代码编写不同，我感觉我的是对的
    return distance


# Function to sort by values
def sort_by_values(list1, values):
    sorted_list = []
    while (len(sorted_list) != len(list1)):
        if index_of(min(values), values) in list1:
            sorted_list.append(index_of(min(values), values))
        values[index_of(min(values), values)] = math.inf
    return sorted_list


# Function to find index of list
def index_of(a, list):
    for i in range(0, len(list)):
        if list[i] == a:
            return i
    return -1


class MA:
    def __init__(self,Pt,Jn,Mn,At,An,Aw,epk,epik,ea_load,ea_unload,ea_idle,pop_size,pc,pm):
        self.fja=FJSP_agv(Pt,Jn,Mn,At,An,Aw,epk,epik,ea_load,ea_unload,ea_idle)
        self.An=An
        self.Pop_size=pop_size
        self.pc=pc
        self.pm=pm
        # self.lr=lr
        self.op_num=[len(Pi) for Pi in self.fja.Pt]
        self.Chromo_list=[]
        for i in range(len (self.op_num)):
            self.Chromo_list.extend([i for _ in range(2*self.op_num[i])])
        self.machine_use = [[] for i in range(len(self.op_num))]
        for j in range(len(self.op_num)):
            for k in range(self.op_num[j]):
                mk = []
                for m in range(self.fja.Mn):
                    if self.fja.Pt[j][k][m] != 9999:
                        mk.append(m)
                self.machine_use[j].append(mk)

    def initial_population(self):
        self.Pop1=[]
        self.Pop2=[]
        self.Pop3=[]
        #工件的初始化规则：
        for i in range(self.Pop_size):
            random.shuffle(self.Chromo_list)
            self.Pop1.append(copy.copy(self.Chromo_list))
            #机器的初始化规则
            machine_order = []
            rand_m=random.uniform(0,1)
            if rand_m<=0.4:#选择规则1，随机初始化
                for j in range(len(self.machine_use)):
                    for k in range(len(self.machine_use[j])):
                        rand_num = random.choice(self.machine_use[j][k])
                        machine_order.append(rand_num)
                self.Pop2.append(machine_order)
            elif 0.4<rand_m<=0.6:#选择规则2，最小的负载
                machine_load = [0 for i in range(self.fja.Mn)]
                for j in range(len(self.machine_use)):
                    for k in range(len(self.machine_use[j])):
                        ma_load = []
                        for m in range(len(self.machine_use[j][k])):
                            ma_load.append(machine_load[self.machine_use[j][k][m]])
                        minimum_value = min(ma_load)
                        indices = []
                        for s in range(len(ma_load)):
                            if ma_load[s] == minimum_value:
                                indices.append(s)
                        logger.info(indices)
                        ma_minindex = random.choice(indices)
                        # ma_minindex=min(range(len(ma_load)), key=lambda i: ma_load[i])
                        machine_load[self.machine_use[j][k][ma_minindex]] += self.fja.Pt[j][k][
                            self.machine_use[j][k][ma_minindex]]
                        machine_order.append(self.machine_use[j][k][ma_minindex])
                self.Pop2.append(machine_order)

            elif 0.6<rand_m<=0.8:#选择规则3，最小的加工时间
                for j in range(len(self.fja.Pt)):
                    for k in range(len(self.fja.Pt[j])):
                        ma_minindex = min(range(len(self.fja.Pt[j][k])),
                                          key=lambda i: self.fja.Pt[j][k][i])  # 找到最小值的索引值
                        machine_order.append(ma_minindex)
                self.Pop2.append(machine_order)
            elif 0.8<rand_m:#选择规则4，最小的完工时间,AGV也需要满足最早运输完成
                AGV_order = [-1 for i in range(sum(self.op_num))]
                machine_order = [-1 for i in range(sum(self.op_num))]
                agv_qs = [0 for i in range(self.fja.Jn)]
                agv_endtime = [0 for i in range(self.An)]
                machine_endtime = [0 for i in range(self.fja.Mn)]
                agv_location = [0 for i in range(self.An)]
                for j in range(len(self.Pop1[i])):
                    agv_qs[self.Pop1[i][j]] += 1
                    if self.Pop1[i][j] == 0:
                        chrom_location = (agv_qs[self.Pop1[i][j]] - 1) // 2
                    else:
                        chrom_location = sum(self.op_num[:self.Pop1[i][j]]) + \
                                         (agv_qs[self.Pop1[i][j]] - 1) // 2
                    if agv_qs[self.Pop1[i][j]] % 2 == 0:  # 表示整除
                        # 表示是去送货,送货跟取货是同一个AGV，由取货来确定agv的选择（工件的完工时间与AGV的结束时间最大值为开始时间）
                        ma_endtime = [0 for i in range(self.fja.Mn)]
                        for m in range(len(machine_endtime)):
                            ma_endtime[m] = max(agv_endtime[AGV_order[chrom_location]] +
                                                self.fja.At[agv_location[AGV_order[chrom_location]]][machine_order[m]] +
                                                self.fja.Pt[self.Pop1[i][j]][(agv_qs[self.Pop1[i][j]] - 1) // 2][m],
                                                machine_endtime[m] +
                                                self.fja.Pt[self.Pop1[i][j]][(agv_qs[self.Pop1[i][j]] - 1) // 2][m])
                        ma_minindex = min(range(len(ma_endtime)), key=lambda i: ma_endtime[i])
                        machine_endtime[ma_minindex] = ma_endtime[ma_minindex]
                        machine_order[chrom_location] = ma_minindex  # 更新机器染色体
                        agv_location[AGV_order[chrom_location]] = ma_minindex  # 更新AGV的位置
                        agv_endtime[AGV_order[chrom_location]] += self.fja.At[agv_location[AGV_order[chrom_location]]][machine_order[m]]  # 更新AGV的完成时间
                    else:
                        a_endtime = [0 for i in range(self.An)]
                        for k in range(len(agv_endtime)):
                            # 表示去取货,确定AGV
                            # 如果第一道工序，则是去LU取货
                            if (agv_qs[self.Pop1[i][j]] == 1):
                                a_endtime[k] = agv_endtime[k] + self.fja.At[agv_location[k]][0]
                            # 否则去前一道工序加工处的位置取货
                            else:  # agv结束上一个任务才可以去取货，取货的结束时间是（AGV到达的时间与工件加工结束时间的最大值）
                                a_endtime[k] = max(
                                    agv_endtime[k] + self.fja.At[agv_location[k]][machine_order[chrom_location]],
                                    machine_endtime[machine_order[chrom_location]])  # AGV到达时间与机器完工时间的最大值#更新AGV的完工时间
                        agv_minindex = min(range(len(a_endtime)), key=lambda i: a_endtime[i])  # 找到最小值的索引值
                        agv_location[agv_minindex] = machine_order[chrom_location]  # 更新AGV的位置
                        agv_endtime[agv_minindex] = a_endtime[agv_minindex]  # 更新AGV的完工时间
                        AGV_order[chrom_location] = agv_minindex  # 更新AGV染色体
                self.Pop2.append(machine_order)
                self.Pop3.append(AGV_order)
                continue#不执行下面的了
            #AGV的初始化规则
            AGV_order =[]
            rand_v = random.uniform(0, 1)
            if rand_v<=0.4:#选择规则1，随机初始化
                AGV_order=[random.randint(0,self.An-1) for i in range(sum(self.op_num))]
                self.Pop3.append(AGV_order)
            elif 0.4<rand_v<=0.6:#选择规则2，最小的任务数，选择最早完工的AGV，要是AGV完工时间一样，则随机选择一个AGV
                AGV_choice=[0 for i in range(self.An)]
                for j in range(sum(self.op_num)):
                    agv_minindex = min(range(len(AGV_choice)), key=lambda i: AGV_choice[i])
                    AGV_order.append(agv_minindex)
                    AGV_choice[agv_minindex]+=1
                self.Pop3.append(AGV_order)
            elif 0.6<rand_v<=0.8:#选择规则3，一个工件选择一个AGV，减少AGV的无效运输
                for s in range(len(self.op_num)):
                    agv_num = random.randint(0, self.An - 1)
                    for j in range(self.op_num[s]):
                        AGV_order.append(agv_num)  # 同一个工件的所有任务都选择一个AGV
                agv_qs = [0 for l in range(self.fja.Jn)]
                agv_location = [0 for l in range(self.An)]
                for j in range(len(self.Pop1[i])):
                    agv_qs[self.Pop1[i][j]] += 1
                    if self.Pop1[i][j] == 0:
                        chrom_location = (agv_qs[self.Pop1[i][j]] - 1) // 2
                    else:
                        chrom_location = sum(self.op_num[:self.Pop1[i][j]]) + (agv_qs[self.Pop1[i][j]] - 1) // 2
                    if agv_qs[self.Pop1[i][j]] % 2 == 0:  # 表示整除，表示去送货
                        agv_location[AGV_order[chrom_location]] = machine_order[chrom_location]  # 更新AGV的位置
                    else:  # 表示去取货，通过取货来确定AGV的选择
                        agv_transtime = [0 for i in range(self.An)]
                        for k in range(len(agv_transtime)):
                            # 表示去取货,确定AGV
                            # 如果第一道工序，则是去LU取货
                            if (agv_qs[self.Pop1[i][j]] == 1):
                                agv_transtime[k] = self.fja.At[agv_location[k]][0]
                            # 否则去前一道工序加工处的位置取货
                            else:  # agv结束上一个任务才可以去取货，取货的结束时间是（AGV到达的时间与工件加工结束时间的最大值）
                                agv_transtime[k] = self.fja.At[agv_location[k]][
                                    machine_order[chrom_location]]  # 更新AGV的完工时间
                        minimum_value = min(agv_transtime)
                        indices = []
                        for s in range(len(agv_transtime)):
                            if agv_transtime[s] == minimum_value:
                                indices.append(s)
                        logger.info(indices)
                        agv_minindex = random.choice(indices)
                        # agv_minindex = min(range(len(agv_transtime)), key=lambda i: agv_transtime[i])  # 找到最小值的索引值
                        agv_location[agv_minindex] = machine_order[chrom_location]  # 更新AGV的位置
                        AGV_order[chrom_location] = agv_minindex  # 更新AGV染色体
                self.Pop3.append(AGV_order)

            elif 0.8<rand_v:#选择规则4，最小的运输时间
                AGV_order = [-1 for i in range(sum(self.op_num))]
                agv_qs = [0 for i in range(self.fja.Jn)]
                agv_location = [0 for i in range(self.An)]
                for j in range(len(self.Pop1[i])):
                    agv_qs[self.Pop1[i][j]] += 1
                    if self.Pop1[i][j] == 0:
                        chrom_location = (agv_qs[self.Pop1[i][j]] - 1) // 2
                    else:
                        chrom_location = sum(self.op_num[:self.Pop1[i][j]]) + (agv_qs[self.Pop1[i][j]] - 1) // 2
                    if agv_qs[self.Pop1[i][j]] % 2 == 0:  # 表示整除，表示去送货
                        agv_location[AGV_order[chrom_location]] = machine_order[chrom_location]  # 更新AGV的位置
                    else:  # 表示去取货，通过取货来确定AGV的选择
                        agv_transtime = [0 for i in range(self.An)]
                        for k in range(len(agv_transtime)):
                            # 表示去取货,确定AGV
                            # 如果第一道工序，则是去LU取货
                            if (agv_qs[self.Pop1[i][j]] == 1):
                                agv_transtime[k] = self.fja.At[agv_location[k]][0]
                            # 否则去前一道工序加工处的位置取货
                            else:  # agv结束上一个任务才可以去取货，取货的结束时间是（AGV到达的时间与工件加工结束时间的最大值）
                                agv_transtime[k] = self.fja.At[agv_location[k]][
                                    machine_order[chrom_location]]  # 更新AGV的完工时间
                        agv_minindex = min(range(len(agv_transtime)), key=lambda i: agv_transtime[i])  # 找到最小值的索引值
                        agv_location[agv_minindex] = machine_order[chrom_location]  # 更新AGV的位置
                        AGV_order[chrom_location] = agv_minindex  # 更新AGV染色体
                self.Pop3.append(AGV_order)
        return self.Pop1,self.Pop2,self.Pop3

    def jbchrom_reset(self, pop1, pop3):  # 调整第一层染色体的顺序，防止AGV的载重超载
        # for Pi in range(len(pop1)):
        c = [0 for i in range(self.fja.Jn)]  # 用来记录装卸及工序号
        c_ = [0 for i in range(self.fja.Jn)]  # 用来记录工序的数量，只有卸载之后才能+1
        a = [0 for i in range(self.fja.An)]  # 用来记录AGV的载重
        a1 = [[] for i in range(self.fja.An)]  # 用来记录AGV的加工工件的路径
        a2 = [[] for i in range(self.fja.An)]  # 用来记录基因的位置
        b = 0  # 需要删除基因的位置
        for Pj in range(len(pop1)):
            j = pop1[Pj]  # 工件号
            if j == 0:
                gx_loca = c_[j]  # 工序号
            else:  # 不是第一个工件
                gx_loca = sum(self.op_num[:j]) + c_[j]  # 工序的位置
            agv = pop3[gx_loca]  # agv号
            if c[j] % 2 == 0:  # 表示该工序是装载，需要判断载重,且为第一道工序
                if a[agv] + 1 > AGV_weight:  # 超过了载重，需要更新基因的位置以及机器的位置
                    # 找到上一个加工的位置的基因 然后调整到目前这个基因的位置前
                    # 先更换位置，然后卸载工序
                    j = a1[agv][-1]  # 上一个工件号
                    for i in range(a2[agv][-1] + 1, len(pop1)):
                        if pop1[i] == j:
                            b = i
                            # logger.info("调整染色体")
                            break
                    pop1.pop(b)  # 删除该位置的索引的值
                    pop1.insert(Pj, j)
                    c_[j] += 1
                    a[agv] -= 1
                else:
                    a[agv] += 1
            else:  # 表示该工序是卸载,需要更新AGV的载重
                c_[j] += 1
                a[agv] -= 1
            c[j] += 1
            a1[agv].append(j)
            a2[agv].append(Pj)
        return pop1

    def cross(self,newpt0):
        new_pt0 = []
        for i in range(self.Pop_size):  # 遍历一半染色体
            t1 = random.randint(0, len(newpt0) - 1)
            t2 = random.randint(0, len(newpt0) - 1)
            while (t1 == t2):
                t2 = random.randint(0, len(newpt0) - 1)
            new_pop1 = [[-1 for n in range(len(newpt0[0][m]))] for m in range(len(newpt0[0]))]
            new_pop2 = [[-1 for n in range(len(newpt0[0][m]))] for m in range(len(newpt0[0]))]
            gj = [i for i in range(self.fja.Jn)]
            if np.random.rand() < self.pc:
                random.shuffle(gj)  # MPX交叉
                # 工序交叉
                split_point = random.randint(1, len(gj) - 1)
                x1 = gj[:split_point]
                x2 = gj[split_point:]
                for j in range(len(newpt0[t1][0])):
                    if newpt0[t1][0][j] in x1:
                        new_pop1[0][j] = newpt0[t1][0][j]
                site = [j for j, var in enumerate(new_pop1[0]) if var == -1]
                e = 0
                for j in range(len(newpt0[t1][0])):
                    if newpt0[t2][0][j] in x2:
                        new_pop1[0][site[e]] = newpt0[t2][0][j]
                        e = e + 1
                for j in range(len(newpt0[t1][0])):
                    if newpt0[t2][0][j] in x2:
                        new_pop2[0][j] = newpt0[t2][0][j]
                site = [j for j, var in enumerate(new_pop2[0]) if var == -1]
                e = 0
                for j in range(len(newpt0[t1][0])):
                    if newpt0[t1][0][j] in x1:
                        new_pop2[0][site[e]] = newpt0[t1][0][j]
                        e = e + 1
                # 机器与AGV的交叉
                random_list = [random.choice([0, 1]) for _ in range(sum(self.op_num))]
                for j in range(len(random_list)):
                    if random_list[j] == 1:
                        new_pop1[1][j] = newpt0[t1][1][j]
                        new_pop1[2][j] = newpt0[t1][2][j]
                        new_pop2[1][j] = newpt0[t2][1][j]
                        new_pop2[2][j] = newpt0[t2][2][j]
                    else:
                        new_pop1[1][j] = newpt0[t2][1][j]
                        new_pop1[2][j] = newpt0[t2][2][j]
                        new_pop2[1][j] = newpt0[t1][1][j]
                        new_pop2[2][j] = newpt0[t1][2][j]
                new_pt0.append(new_pop1)
                new_pt0.append(new_pop2)
            else:
                new_pt0.append(newpt0[t1])
                new_pt0.append(newpt0[t2])
        return new_pt0

    def flatten_3d_to_2d(self,lst):
        flattened_list = []
        for sublist in lst:
            for item in sublist:
                flattened_list.append(item)
        return flattened_list

    def mutation(self, pop):
        new_pop = copy.copy(pop)
        for i in range(len(pop)):
            if np.random.rand() < self.pm:
                random_positions = random.sample(range(len(new_pop[i][0])), 2)
                b = new_pop[i][0][random_positions[0]]
                new_pop[i][0][random_positions[0]] = new_pop[i][0][random_positions[1]]
                new_pop[i][0][random_positions[1]] = b
                random_positions1 = random.sample(range(len(new_pop[i][2])), 3)
                for j in random_positions1:
                    nested_list_3d = copy.copy(self.machine_use)
                    flattened_list_2d = ma.flatten_3d_to_2d(nested_list_3d)
                    rand_num = random.choice(flattened_list_2d[j])
                    new_pop[i][1][j] = rand_num
                    new_pop[i][2][j] = random.choice([s for s in range(self.An) if s != new_pop[i][2][j]])
        return new_pop
    def fitness(self,pop):
        value=[]
        fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number, AGV_weight, epk,
                       epik, ea_load, ea_unload, ea_idle)
        self.Jobs, self.Machines, self.AGVs = fja.decode(pop[0], pop[1], pop[2])

        v1, v2 = fja.total_value()
        value.append(v1)
        value.append(v2)
        return value
    def non_dominate(self,pop):
        # 计算适应度值
        func_value1 = []
        func_value2 = []
        for i in range(len(pop)):
            value = ma.fitness(pop[i])
            func_value1.append(value[0])
            func_value2.append(value[1])
        non_dominated_sorted_solution = fast_non_dominated_sort(func_value1, func_value2)  # 非支配排序
        return non_dominated_sorted_solution, func_value1, func_value2

    def VNS(self,pt):
        fja = FJSP_agv( P_time, Job_number, Machine_number, AGV_Trans,AGV_number,AGV_weight,epk,epik,ea_load,ea_unload,ea_idle)
        self.Jobs, self.Machines, self.AGVs = fja.decode(pt[0], pt[1], pt[2])
        v1, v2 = fja.total_value()
        cri_path,cri_op=fja.critical_path(v1)
        cri_op1=list(chain.from_iterable(cri_op))#将关键工序展开成一维列表
        rule=random.randint(1,5)#随机选择一个规则
        if rule==1:#合并非关键工序的取货，在AGV位置上，减少能耗，这个没有实现
            cri_block = [sublist for sublist in cri_op if len(sublist) > 1]
            if cri_block == []:  # 不存在关键块
                pass
            else:
                random_number = random.choice(cri_block)
                if (random_number[0] // 10 == random_number[1] // 10):  # 如果两个工序是同一个工件的
                    pass
                else:
                    indices0 = [i for i, val in enumerate(pt[0]) if
                                val == random_number[0] // 10 - 1]  # 找到了第一个工序的索引位置
                    indices1 = [i for i, val in enumerate(pt[0]) if
                                val == random_number[1] // 10 - 1]  # 找到了第二个工序的索引位置
                    pt[0][indices0[2 * random_number[0] % 10 - 1]] = random_number[1] // 10 - 1
                    pt[0][indices1[2 * random_number[1] % 10 - 1]] = random_number[
                                                                           0] // 10 - 1  # 交换两个位置的元素。保证第二个工序的放货在第一个的前面
        elif rule==2:#随机选择一个关键的工序，将该工序的放货顺序随机往前插入，更改工序染色体
            random_number = random.choice(cri_op1)
            gj_num = random_number // 10
            gx_num = random_number % 10
            indices=[i for i,val in enumerate(pt[0]) if val==gj_num-1]#找到了工序的索引位置
            if (indices[2*gx_num-1]==indices[2*gx_num-2]+1):#如果两个位置是相挨着的
                pass
            else:
                insert_index=random.randint(indices[2*gx_num-2]+1,indices[2*gx_num-1]-1)
                pt[0].pop(indices[2*gx_num-1])#移除最后一个位置
                pt[0].insert(insert_index,gj_num-1)

        elif rule==3:#随机选择一个关键工序，更换该工序的机器，更改机器染色体
            random_number = random.choice(cri_op1)
            gj_num = random_number // 10
            gx_num = random_number % 10
            change_ma = sum(self.op_num[:gj_num - 1]) + gx_num - 1  # 染色体上的位置
            # 更改染色体上的machine
            machine_num=pt[1][change_ma]
            try:
                candidates=[i for i in range(self.fja.Mn) if self.fja.Pt[gj_num-1][gx_num-1][i]!=9999]#为什么会出错 不是有个if 啥的 error 可以打印出来吗
            except Exception as err:
                logger.info(err)
                traceback.logger.info_exc()  # 打印堆栈信息
            pt[1][change_ma]=random.choice(candidates)
        elif rule==4:#找到关键块，交换头部与尾部的顺序，更换工序染色体
            cri_block=[sublist for sublist in cri_op if len(sublist)>1]
            if cri_block==[]:#不存在关键块
                pass
            else:
                random_number=random.choice(cri_block)
                if (random_number[0]//10==random_number[1]//10):#如果两个工序是同一个工件的
                    pass
                else:
                    indices0=[i for i,val in enumerate(pt[0]) if val==random_number[0]//10-1]#找到了第一个工序的索引位置
                    indices1=[i for i,val in enumerate(pt[0]) if val==random_number[1]//10-1]#找到了第二个工序的索引位置
                    pt[0][indices0[2*random_number[0]%10-1]]=random_number[1]//10-1
                    pt[0][indices1[2*random_number[1]%10-1]]=random_number[0]//10-1#交换两个位置的元素。保证第二个工序的放货在第一个的前面

        else:#随机选择一个关键工序，更换AGV
            random_number=random.choice(cri_op1)
            gj_num=random_number//10
            gx_num=random_number%10
            change_agv=sum(self.op_num[:gj_num-1])+gx_num-1#染色体上的位置
            #更改染色体上的AGV
            agv_num=pt[2][change_agv]
            candidates=[i for i in range(self.An) if i!=agv_num]
            pt[2][change_agv]=random.choice(candidates)
        return pt

    def gantt(self,pop1,pop2,pop3):
        fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number, AGV_weight, epk, epik, ea_load,
                       ea_unload, ea_idle)
        self.Jobs, self.Machines, self.AGVs = fja.decode(pop1, pop2, pop3)
        fja.Gatt()

    def sav_energy(self,pop1,pop2,pop3):#节约能耗，甘特图整体右移，时间不变,再求一遍能耗目标值
        #关键工序不可以动，只能移动非关键工序
        fja = FJSP_agv( P_time, Job_number, Machine_number, AGV_Trans,AGV_number,AGV_weight,epk,epik,ea_load,ea_unload,ea_idle)
        self.Jobs, self.Machines, self.AGVs = fja.decode(pop1, pop2, pop3)
        v1, v2 = fja.total_value()
        cri_path, cri_op = fja.critical_path(v1)
        cri_op1 = list(chain.from_iterable(cri_op))  # 将关键工序展开成一维列表
        last_op_m=[]#每台机器上的最后一道工序
        for j in range(self.fja.Mn):
            Machine =self.Machines[j]
            if(len(Machine.M_end)>0):
                last_op_m.append(Machine.M_job[-1])
        gx_num=[]
        for i in range(len(self.op_num)):
            for j in range(self.op_num[i]):
                gx_num.append((i+1)*10+(j+1))
        gx_num1=copy.copy(gx_num)
        for i in gx_num:
            if i in cri_op1 or i in last_op_m:#如果是关键工序，则不能改变时间,或者是机器上的最后一道工序不动，只运动其他工序
                gx_num1.remove(i)
        while(len(gx_num1)>0):
            random_gj=random.choice(gx_num1)
            gx_num1.remove(random_gj)
            random_gj_num = random_gj // 10#工件号
            random_gx_num = random_gj % 10#工序号
            ma_num=pop2[sum(self.op_num[:random_gj_num-1])+random_gx_num-1]
            Machine=self.Machines[ma_num]
            index=[i for i,var in enumerate (Machine.M_job) if var==random_gj]
            index0=index[0]
            gj1=Machine.M_job[index0+1]#同一个机器上的后边一个工件
            if (random_gj+1) in gx_num:#不是工件的最后一道工序，判定下一道工序是否时间确定
                if (random_gj+1) not in gx_num1 and gj1 not in gx_num1:#工件的后一道工序已经确定时间，同一机器上的后一道工序确定时间
                    #当前工序的下一道工序的时间确定了
                    agv_num_last=pop3[sum(self.op_num[:random_gj_num-1])+random_gx_num]
                    agv=self.AGVs[agv_num_last]
                    index1=[i for i,val in enumerate(agv.A_job) if val==random_gj+1]
                    if index1==[]:
                        pass
                    else:
                        if(Machine.M_end[index0]==Machine.M_start[index0+1]) :#这边有问题，不知道原因
                            pass
                        elif(Machine.M_end[index0]==agv.A_start[index1[1]]):
                            pass
                        else:
                            pro = Machine.M_end[index0] - Machine.M_start[index0]
                            Machine.M_end[index0] = min(Machine.M_start[index0 + 1],agv.A_start[index1[1]])#取后边两者的最小值
                            Machine.M_start[index0] = Machine.M_end[index0] - pro  # 更新了工序的开始时间
                else:
                    pass

            else:#是工件的最后一道工序，只需考虑该工件在机器后边的那道工序的时间是否确定，若没有，pass,确定了就可以确定当前工序的时间
                if gj1 in gx_num1:#这个工序的时间还没有确定
                    pass
                else:#机器上后边的工件的时间已经确定
                    if(Machine.M_end[index0]==Machine.M_start[index0+1]):
                        pass#没有空隙
                    else:
                        pro=Machine.M_end[index0]-Machine.M_start[index0]
                        Machine.M_end[index0]=Machine.M_start[index0+1]
                        Machine.M_start[index0]=Machine.M_end[index0]-pro#更新了工序的开始时间
        #更新完成之后再去求一次能耗目标值
        v1, v2 = fja.total_value()
        return v1,v2
    def remove_duplicates(self,lists):
        unique_lists = []
        duplicates = set()

        for i in range(len(lists)):
            if i not in duplicates:
                unique_lists.append(lists[i])
                for j in range(i + 1, len(lists)):
                    if lists[i] == lists[j]:
                        duplicates.add(j)
        return unique_lists

    def Mating(self, pop_num, rank,crowd):#锦标赛选择
        Matingpool=[]
        num_range = list(range(0, pop_num))
        for b in range(pop_num):
            #随机选择两个个体，进行比较，好的个体留下,避免重复选择两个个体
            select_numbers=random.sample(num_range,2)
            p=select_numbers[0]
            q=select_numbers[1]
            if rank[p] < rank[q]:
                best = p
            elif rank[p] > rank[q]:
                best = q
            elif crowd[p] > crowd[q]:
                best = p
            else:
                best = q
            Matingpool.append(best)
        return Matingpool
    def Pareto_crowd(self,non_dominated_sorted_solution,func1,func2):
        crowding_distance_values = []
        for j in range(0, len(non_dominated_sorted_solution)):
            crowding_distance_values.append(
                crowding_distance(func1[:], func2[:], non_dominated_sorted_solution[j][:]))  # 拥挤度计算
        # 存储层级rank与拥挤度crowd
        rank = [-1 for i in range(len(func1))]
        crowd = [-1 for i in range(len(func1))]
        for j in range(0, len(non_dominated_sorted_solution)):
            for k in range(0, len(non_dominated_sorted_solution[j])):
                rank[non_dominated_sorted_solution[j][k]] = j
                crowd[non_dominated_sorted_solution[j][k]] = crowding_distance_values[j][k]
        new_solution = []
        for i in range(0, len(non_dominated_sorted_solution)):
            non_dominated_sorted_solution2_1 = [
                index_of(non_dominated_sorted_solution[i][j], non_dominated_sorted_solution[i]) for j in
                range(0, len(non_dominated_sorted_solution[i]))]  # 返回每一层的最小值的索引，这个写的太繁琐，其实就是变成，0，1，2，3索引
            front22 = sorted(non_dominated_sorted_solution2_1, key=lambda x: crowding_distance_values[i][x])  # 根据拥挤度距离的大小进行排序，先小后打，索引需要反转
            front = [non_dominated_sorted_solution[i][front22[j]] for j in
                     range(0, len(non_dominated_sorted_solution[i]))]
            front.reverse()
            for value in front:
                new_solution.append(value)
                if (len(new_solution) == self.Pop_size):  # 其实就是根据排序和拥挤度取到足够的个体
                    break
            if (len(new_solution) == self.Pop_size):
                break
        return new_solution
    def dominated_by(self, func_value1, func_value2):
        if (func_value1[0]<=func_value2[0] and func_value1[1]<=func_value2[1]):
            return 0
        elif(func_value2[0]<=func_value1[0] and func_value2[1]<=func_value1[1]):
            return 1
        else:
            return -1
    def main(self,pt0):
        # 求适应度
        # func_value1 = []
        # func_value2 = []
        # logger.info(f'解码')
        # for i in range(len(pop1)):
        #     fja = FJSP_agv( P_time, Job_number, Machine_number, AGV_Trans,AGV_number,AGV_weight,epk,epik,ea_load,ea_unload,ea_idle)
        #     self.Jobs, self.Machines, self.AGVs = fja.decode(pt0[i][0], pt0[i][1], pt0[i][2])
        #     v1, v2 = fja.total_value()
        #     func_value1.append(v1)
        #     func_value2.append(v2)
        # logger.info(f'非支配排序')
        # non_dominated_sorted_solution = fast_non_dominated_sort(func_value1, func_value2)  # 非支配排序
        # crowding_distance_values = []
        # logger.info(f'拥挤度计算')
        # for j in range(0, len(non_dominated_sorted_solution)):
        #     crowding_distance_values.append(
        #         crowding_distance(func_value1[:], func_value2[:], non_dominated_sorted_solution[j][:]))  # 拥挤度计算
        # logger.info(f'拥挤度计算end')
        # # 存储层级rank与拥挤度crowd
        # rank = [-1 for i in range(len(func_value1))]
        # crowd = [-1 for i in range(len(func_value1))]
        # for j in range(0, len(non_dominated_sorted_solution)):
        #     for k in range(0, len(non_dominated_sorted_solution[j])):
        #         rank[non_dominated_sorted_solution[j][k]] = j
        #         crowd[non_dominated_sorted_solution[j][k]] = crowding_distance_values[j][k]
        # # 锦标赛选择个体
        # logger.info(f'锦标赛选择个体')
        # new_solution = ma.Mating(len(func_value1), rank, crowd)
        # new1_pt0=[]
        # for n in new_solution:
        #     new1_pt0.append(pt0[n])
        # logger.info(f'交叉')
        new2_pt0 = ma.cross(pt0)#交叉
        logger.info(f'修整工序染色体')
        for i in range(len(new2_pt0)):
            pt0_1 = ma.jbchrom_reset(new2_pt0[i][0], new2_pt0[i][2])
            new2_pt0[i][0] = pt0_1 # 修整工序染色体
        logger.info(f'变异')  
        new3_pt0 = ma.mutation(pt0)#变异
        logger.info(f'修整工序染色体2')
        for i in range(len(new3_pt0)):
            pt0_1 = ma.jbchrom_reset(new3_pt0[i][0], new3_pt0[i][2])
            new3_pt0[i][0] = pt0_1 # 修整工序染色体  # 修整工序染色体
        logger.info(f'修整工序染色体2 end')  
        # 交叉变异
        # 修整染色体
        logger.info(f'选择1/4个体进行VNS')  
        #选择1/4个体进行VNS
        pt = pt0+ new2_pt0 + new3_pt0
        pt = ma.remove_duplicates(pt)
        new_list,func1,func2=ma.non_dominate(pt)
        R1 = []
        R2 = []
        new_solution = ma.Pareto_crowd(new_list, func1, func2)
        for n in new_solution:
            R1.append(pt[n])
        #######

        for p in range(1):
            for j in new_list[p]:
                # 前沿个体采用关键路径规则进行优化
                POP = ma.VNS(pt[j])
                pt0_1 = ma.jbchrom_reset(POP[0], POP[2])  # 修复染色体
                POP[0] = pt0_1
                R2.append(POP)


        R3 = R1 + R2
        pt0 = []
        logger.info(f'non_dominate')
        non_dominated_sorted_solutio1, func11, func21 = ma.non_dominate(R3)
        new_solution = ma.Pareto_crowd(non_dominated_sorted_solutio1, func11, func21)
        for n in new_solution:
            pt0.append(R3[n])
        logger.info(f'non_dominate end')
        objective_1=[func11[i] for i in new_solution]
        objective_2 = [func21[i] for i in new_solution]
        return pt0,objective_1,objective_2



if __name__=="__main__":
    prefix = " "
    try:
        prefix = sys.argv[1] + "_"
        # logger.info("prefix = %s", prefix)
    except Exception as e:
        pass
        # logger.warning("没有系统参数,建议直接修改prefix")

    # ps=[50,10,150,200]
    # pc=[0.6,0.7,0.8,0.9]
    # pm=[0.05,0.1,0.15,0.2]
    # Maxit=[10,200,300,400]
    # parameter=[[1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4],
    #            [1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4],
    #            [1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1],
    #            [1,2,3,4,3,4,1,2,4,3,2,1,2,1,4,3],
    #            [1,2,3,4,4,3,2,1,2,1,4,3,3,4,1,2]]
    # for t in range(16):
    #     logger.info('参数次数******', t)
    #     ps1=ps[parameter[0][t]-1]
    #     pc1=pc[parameter[1][t]-1]
    #     pm1=pm[parameter[2][t]-1]
    #     Maxit1=Maxit[parameter[3][t]-1]
    #     # lr1=lr[parameter[4][t]-1]
    # ps1=200
    # pc1=0.8
    # pm1=0.2
    # Maxit1=400
    # lr1=0.2
    # files=[1,3,7,9,11,13,15]
    # for frequence in range(2):
        for file_index in range(9,10):
            # logger.info(f'../1_Brandimarte/BrandimarteMk{file_index}.fjs')
            # P_time, Job_number, Machine_number, AGV_Trans, AGV_number,AGV_weight, epk, epik, ea_load, ea_unload, ea_idle = Get_Problem(
            #     f'../1_Brandimarte/BrandimarteMk{file_index}.fjs')

            # P_time, Job_number, Machine_number, AGV_Trans, AGV_number, AGV_weight, epk, epik, ea_load, ea_unload, ea_idle = Get_Problem(
            #     r'C:\ahua\博士文件夹\算法代码\GA-AGV-uncompetition\1_Brandimarte\BrandimarteMk1.fjs')
            GEN=400
            ps1=200
            pc1=0.8
            pm1=0.2

            Machine_number=10
            Job_number=7
            AGV_number = 3
            AGV_weight = 2
            # epk = [random.randint(10, 18) for i in range(Machine_number)]
            epk=[10,8,12,10,7,13,11,8,10,11]
            # epik = [random.randint(1, 4) for i in range(Machine_number)]
            epik=[2,1,2,1,2,1,1,2,3,2]
            ea_unload = 3.5
            ea_load = 0.5
            ea_idle = 0.5
            P_time = [[
                [9999, 9999, 6, 4, 9999, 9999, 9999, 9999, 9999, 9999],
                [2, 3, 9999, 9999, 4, 9999, 9999, 9999, 9999, 9999],
                [4, 3, 3, 6, 9999, 4, 9999, 9999, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 4, 3, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 3, 5]
            ], [
                [9999, 9999, 3, 4, 9999, 9999, 9999, 9999, 9999, 9999],
                [2, 2, 4, 3, 9999, 9999, 9999, 9999, 9999, 9999],
                [3, 2, 9999, 9999, 4, 9999, 9999, 9999, 9999, 9999],
                [9999,   9999, 4, 6, 9999, 9999, 9999, 9999, 9999, 9999],
                [3, 3, 9999, 9999, 9999, 9999, 9999, 9999, 2, 4],
                [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 7, 8]

            ], [
                [5, 7, 9999, 9999, 9999, 9999, 9999, 9999, 9, 6],
                [3, 5, 9999, 9999, 3, 9999, 9999, 9999, 9999, 9999],
                [7, 5, 9999, 9999, 9999, 4, 9999, 9999, 9999, 9999],
                [3, 2, 2, 4, 9999, 9999, 9999, 9999, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 6, 9, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 4, 9999, 9999, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 4, 7, 9999, 9999]
            ], [
                [4, 4, 3, 6, 9999, 9999, 9999, 9999, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9, 7],
                [2, 1, 9999, 9999, 3, 9999, 9999, 9999, 9999, 9999],
                [9999, 9999, 7, 6,7, 9999, 9999, 9999, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 4, 9999, 9999, 9999, 9999],
                [3, 3, 9999, 9999, 9999, 9999, 9999, 9999, 4, 5]
            ], [
                [9999, 9999, 5, 6, 5, 9999, 9999, 9999, 9999, 9999],
                [3, 2, 9999, 9999, 2, 9999, 9999, 9999, 9999, 9999],
                [6, 6, 5, 8, 5, 9999, 9999, 9999, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 2, 3],
                [9999, 9999, 9999, 9999, 9999, 9999, 4, 3, 9999, 9999],
                [3, 2, 9999, 9999, 2, 9999, 9999, 9999, 9999, 9999]

            ], [
                [3, 4, 9999, 9999, 9999, 9999, 9999, 9999, 5, 5],
                [9999, 9999, 9999, 9999, 9999, 9999, 4, 3, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 3, 9999, 9999, 9999, 9999],
                [9999, 9999, 3, 4, 3, 9999, 9999, 9999, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 4, 6],
                [3, 3, 2, 4, 9999, 9999, 9999, 9999, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 3, 5, 9999, 9999]
            ], [
                [3, 3, 9999, 2, 9999, 9999, 9999, 9999, 4, 5],
                [3, 2, 9999, 9999, 3, 9999, 9999, 9999, 9999, 9999],
                [1, 9999, 2, 1, 9999, 9999, 9999, 9999, 9999, 9999],
                [9999, 9999, 4, 4, 3, 9999, 9999, 9999, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 3, 4, 9999, 9999],
                [9999, 9999, 9999, 9999, 9999, 9999, 5, 4, 9999, 9999],
                [3, 4, 9999, 4, 9999, 9999, 9999, 9999, 3, 5]
            ]]

            AGV_Trans = [[0, 1.5, 3.9, 3.1, 3.3, 4.1, 4.6, 3.3, 3.7, 3.0, 3.3],
                         [1.5, 0, 2.4, 1.6, 1.8, 2.4, 3.1, 1.8, 2.2, 1.5, 1.8],
                         [3.9, 2.4, 0, 1.7, 1.5, 2.3, 2.2, 1.2, 1.0, 2.5, 2.2],
                         [3.1, 1.6, 1.7, 0, 1.0, 2.2, 2.8, 1.9, 1.2, 1.4, 1.2],
                         [3.3, 1.8, 1.5, 1.0, 0, 2.7, 2.4, 1.8, 1.6, 2.0, 1.8],
                         [4.1, 2.4, 2.3, 2.2, 2.7, 0, 1.8, 2.8, 2.6, 1.2, 3.0],
                         [4.6, 3.1, 2.2, 2.8, 2.4, 1.8, 0, 2.0, 2.3, 2.8, 2.6],
                         [3.3, 1.8, 1.2, 1.9, 1.8, 2.8, 2.0, 0, 1.0, 1.2, 1.4],
                         [3.7, 2.2, 1.0, 1.2, 1.6, 2.6, 2.3, 1.0, 0, 1.5, 1.7],
                         [3.0, 1.5, 2.5, 1.4, 2.0, 1.2, 2.3, 1.2, 1.5, 0, 1.1],
                         [3.3, 1.8, 2.2, 1.2, 1.8, 3.0, 2.6, 1.4, 1.7, 1.1, 0]
                         ]

            # for j1 in range(len(AGV_W)):
            #     logger.info(f'AGV的载重是多少************ {AGV_W[j1]}')
            # AGV_weight=AGV_W[2]
            # for i1 in range(len(AGV_N)):
            #     logger.info(f'AGV的数量是多少************ {AGV_N[i1]}')
            #     AGV_number = AGV_N[i1]

            ma = MA( P_time, Job_number, Machine_number, AGV_Trans,AGV_number,AGV_weight,epk,epik,ea_load,ea_unload,ea_idle,ps1,pc1,pm1)
            pop1, pop2, pop3 = ma.initial_population()  # 初始化
            for j in range(len(pop1)):
                pop = ma.jbchrom_reset(pop1[j], pop3[j])  # 修整工序染色体
                pop1[j]=pop
            pt0 = []
            for i in range(len(pop1)):
                pt1 = []
                pt1.append(pop1[i])
                pt1.append(pop2[i])
                pt1.append(pop3[i])
                pt0.append(pt1)
            pareto_front_solution=[]
            func1_best=[]
            func2_best=[]
            best_pop=[]
            for g in range(GEN):
                logger.info(f'迭代次数******{g}')
                pt0,func_value1,func_value2=ma.main(pt0)
                logger.info(f'初始化种群')
                #更新迭代完成后，进行能耗的优化
                # func1_object = func_value1
                # func2_object = func_value2
                func1_object=[]
                func2_object=[]
                for j in range(len(pt0)):
                    value1,value2=ma.sav_energy( pt0[j][0],pt0[j][1],pt0[j][2])
                    if (value1<=func_value1[j]) and(value2<=func_value2[j]):
                        func1_object.append(value1)#这边的完工时间不应该会变化的
                        func2_object.append(value2)
                    else:
                        func1_object.append(func_value1[j])
                        func2_object.append(func_value2[j])
                #进行非支配排序，找到非支配解，并把非支配解的目标保存
                logger.info(f'非支配排序')
                non_dominated_sorted_solution = fast_non_dominated_sort(func1_object, func2_object)  # 非支配排序
                if g==GEN-1:
                    for j in non_dominated_sorted_solution[0]:
                        print(pt0[j])
                        print(func1_object[j])
                        print(func2_object[j])
                        ma.gantt(pt0[j][0],pt0[j][1],pt0[j][2])
                objectives=[]
                for j in non_dominated_sorted_solution[0]:
                    object=[]
                    object.append(func1_object[j])
                    object.append(func2_object[j])
                    objectives.append(object)
                pareto_front_solution.append(objectives)
                ob=np.asarray(objectives)
                func1_best.append(min(ob[:,0]))
                func2_best.append(min(ob[:,1]))

                print(func1_best)
                print(func2_best)
                # logger.info(pareto_front_solution)
                # logger.info(func1_best)

                logger.info(func2_best)
            save_plot_main_data("pic_data/FJSP_MlAGV/data/", pareto_front_solution,func1_best,func2_best,GEN,prefix=prefix)

    #画图

#导入数据案例
#导出数据，保存，画出与其他算法的迭代对比图，三种
#计算IGD,SC，策略验证实验IGD,SC


#参数验证实验,popsize,pc,pm,ls,迭代次数
#MILP验证实验
#AGV的数量实验，多载的数量1，2，3，4


#伪代码两个







    #非支配排序得到前1/4个体进行VNS
    #VNS
    #修整染色体
    #求解目标值，完工时间与能耗
    #优化能耗，求解目标值，完工时间与能耗

    # p1,p2,p3=ma.initial_population()#初始化
    # jb=ma.jbchrom_reset()#调整染色体顺序

    # '''logger.info(P_time)
    # logger.info(ma.op_num)
    # logger.info(ma.machine_use)
    # logger.info(ma.initial_population())
    # # 这个函数有几个返回值  就可以用几个变量接收返回值
    # p1,p2,p3 = ma.initial_population()
    # logger.info("我是p1:",p1)
    # logger.info("我是p2:", p2)
    # logger.info("我是p3:", p3)
    # ### 这里还有一个概念  就是函数返回值是个tuple
    # tup =  ma.initial_population()  #
    # logger.info("我是p3,通过这个tuple加索引的方式调用：", tup[2])
    #
    # #  因为 pop1 pop2 pop3 都是 class MA 中的 self 调用的  就是 MA 的实例 ma 的属性值   明白了吗？ok
    # logger.info("我是pop3：通过实例ma打印：",ma.Pop3)'''



