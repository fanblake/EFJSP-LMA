#XU 多目标
"""Multi-objective green scheduling of integrated flexible job shop and
automated guided vehicles"""
"""
1、初始化种群P（大小为N）和空的外部储备集P'（大小为N'）
2、将P中的非支配解复制到P'
3、如果P'大小超过N'，通过聚类的方法移除
    1、初始化聚类集合C，P'中每个个体构成一个聚类
    2、如果C《N'（即簇的所有个体小于给定的最大值）去5，否则，3
    3、计算不同簇间的两个个体的距离
    4、将存在最小距离的两个簇合并，去2
    5、从每个簇中选择一个具有代表的个体。（簇内某个个体到其他个体的平均距离最小的个体）作为代表解
4、计算P和P'中所有个体的适应值
5、从P+P'中进行选择操作，直到交配池填满
6、交配池中个体进行变异和交叉操作
7、是否达到最大迭代次数，返回第二步
"""
from Data_transfer import Get_Problem
from FJSP_agv import FJSP_agv,Job,Machine,Agv
import random
import copy
import numpy as np
from itertools import chain
from myplot import save_plot_main_data
import math
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
                max(values2) - min(values2))  # 这边我进行了修整，不知道对不对，跟参考的代码编写不同，我感觉我的是对的
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
class EHA:
    def __init__(self,Pt,Jn,Mn,At,An,Aw,epk,epik,ea_load,ea_unload,ea_idle,pop_size=200,pc=0.8,pm=0.2):
        self.fja=FJSP_agv(Pt,Jn,Mn,At,An,Aw,epk,epik,ea_load,ea_unload,ea_idle)
        self.An=An
        self.Pop_size=pop_size
        self.pc=pc
        self.pm=pm
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
        self.archive_size=50

    def initial_population(self):  # 随机初始化规则,文章采用了全局，局部，随机初始化规则，我只采用了随机初始化规则
        Pop1 = []
        Pop2 = []
        Pop3 = []
        # 工件的初始化规则：
        for i in range(self.Pop_size):
            random.shuffle(self.Chromo_list)
            Pop1.append(copy.copy(self.Chromo_list))
            machine_order = []
            for j in range(len(self.machine_use)):
                for k in range(len(self.machine_use[j])):
                    rand_num = random.choice(self.machine_use[j][k])
                    machine_order.append(rand_num)
            Pop2.append(machine_order)
            AGV_order = [random.randint(0, self.An - 1) for i in range(sum(self.op_num))]
            Pop3.append(AGV_order)
        return Pop1, Pop2, Pop3

    def cross(self,newpt0):
        new_pt0 = []
        for i in range(self.Pop_size // 2):  # 遍历一半染色体
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
                    flattened_list_2d = EH.flatten_3d_to_2d(nested_list_3d)
                    rand_num = random.choice(flattened_list_2d[j])
                    new_pop[i][1][j] = rand_num
                    new_pop[i][2][j] = random.choice([s for s in range(self.An) if s != new_pop[i][2][j]])
        return new_pop

    def jbchrom_reset(self,pop1,pop3):#调整第一层染色体的顺序，防止AGV的载重超载
        # for Pi in range(len(pop1)):
        c = [0 for i in range(self.fja.Jn)]  # 用来记录装卸及工序号
        c_ = [0 for i in range(self.fja.Jn)]  # 用来记录工序的数量，只有卸载之后才能+1
        a = [0 for i in range(self.fja.An)]  # 用来记录AGV的载重
        a1 = [[] for i in range(self.fja.An)]  # 用来记录AGV的加工工件的路径
        a2 = [[] for i in range(self.fja.An)]  # 用来记录基因的位置
        b = 0  # 需要删除基因的位置
        for Pj in range(len(pop1)):
            j=pop1[Pj]#工件号
            if j==0:
                gx_loca=c_[j]#工序号
            else:#不是第一个工件
                gx_loca=sum(self.op_num[:j])+c_[j]#工序的位置
            agv=pop3[gx_loca]#agv号
            if c[j]%2==0: #表示该工序是装载，需要判断载重,且为第一道工序
                if a[agv]+1>AGV_weight:#超过了载重，需要更新基因的位置以及机器的位置
                    #找到上一个加工的位置的基因 然后调整到目前这个基因的位置前
                    #先更换位置，然后卸载工序
                    j=a1[agv][-1]#上一个工件号
                    for i in range(a2[agv][-1]+1,len(pop1)):
                        if pop1[i]==j:
                            b=i
                            # print("调整染色体")
                            break
                    pop1.pop(b)#删除该位置的索引的值
                    pop1.insert(Pj,j)
                    c_[j] += 1
                    a[agv] -= 1
                else:
                    a[agv]+=1
            else:#表示该工序是卸载,需要更新AGV的载重
                c_[j]+=1
                a[agv]-=1
            c[j] += 1
            a1[agv].append(j)
            a2[agv].append(Pj)
        return pop1
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
            value=EH.fitness(pop[i])
            func_value1.append(value[0])
            func_value2.append(value[1])
        non_dominated_sorted_solution = fast_non_dominated_sort(func_value1, func_value2)  # 非支配排序
        return non_dominated_sorted_solution,func_value1,func_value2
    def VNS(self,POP):
        fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk,
                       epik, ea_load, ea_unload, ea_idle)
        self.Jobs, self.Machines, self.AGVs = fja.decode(POP[0], POP[1], POP[2])
        v1, v2 = fja.total_value()
        cri_path,cri_op=fja.critical_path(v1)
        cri_op1=list(chain.from_iterable(cri_op))#将关键工序展开成一维列表
        rule=random.randint(1,2)#随机选择一个规则
        if rule==1:#找到关键块，交换头部与尾部的顺序，更换工序染色体
            cri_block=[sublist for sublist in cri_op if len(sublist)>1]
            if cri_block==[]:#不存在关键块
                pass
            else:
                random_number=random.choice(cri_block)
                if (random_number[0]//10==random_number[1]//10):#如果两个工序是同一个工件的
                    pass
                else:
                    indices0=[i for i,val in enumerate(POP[0]) if val==random_number[0]//10-1]#找到了第一个工序的索引位置
                    indices1=[i for i,val in enumerate(POP[0]) if val==random_number[1]//10-1]#找到了第二个工序的索引位置
                    POP[0][indices0[2*random_number[0]%10-1]]=random_number[1]//10-1
                    POP[0][indices1[2*random_number[1]%10-1]]=random_number[0]//10-1#交换两个位置的元素。保证第二个工序的放货在第一个的前面

        else:#随机选择一个关键工序，更换该工序的机器，更改机器染色体
            random_number = random.choice(cri_op1)
            gj_num = random_number // 10
            gx_num = random_number % 10
            change_ma = sum(self.op_num[:gj_num - 1]) + gx_num - 1  # 染色体上的位置
            # 更改染色体上的machine
            machine_num=POP[1][change_ma]
            candidates=[i for i in range(self.fja.Mn) if self.fja.Pt[gj_num-1][gx_num-1][i]!=9999]
            POP[1][change_ma]=random.choice(candidates)
        return POP

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
    def main(self,pt0):
        pareto_front_solution = []
        func1_best = []
        func2_best = []
        for g in range(GEN):
            print("迭代次数：",g)
            new_pt0 = EH.cross(pt0)#交叉
            for i in range(len(new_pt0)):
                pt0_1 = EH.jbchrom_reset(new_pt0[i][0], new_pt0[i][2])
                new_pt0[i][0] = pt0_1
            new_pt1 = EH.mutation(pt0)#变异
            for i in range(len(pt0)):
                pt0_1 = EH.jbchrom_reset(new_pt1[i][0], new_pt1[i][2])
                new_pt1[i][0] = pt0_1
            pt1=pt0+new_pt0+new_pt1
            non_dominated_sorted_solution,func1,func2=EH.non_dominate(pt1)
            R1=[]
            R2=[]
            new_solution=EH.Pareto_crowd(non_dominated_sorted_solution,func1,func2)
            for n in new_solution:
                R1.append(pt1[n])
            #######
            for j in non_dominated_sorted_solution[0]:
                #前沿个体采用关键路径规则进行优化
                POP = EH.VNS(pt1[j])
                pt0_1 = EH.jbchrom_reset(POP[0], POP[2])#修复染色体
                POP[0] = pt0_1
                R2.append(POP)
            R3=R1+R2
            pt0=[]
            non_dominated_sorted_solutio1,func11,func21 = EH.non_dominate(R3)
            new_solution = EH.Pareto_crowd(non_dominated_sorted_solutio1,func11,func21)
            for n in new_solution:
                pt0.append(R3[n])
            #返回前沿个体的所有解
            objectives=[]
            for i in non_dominated_sorted_solutio1[0]:
                object = []
                object.append(func11[i])
                object.append(func21[i])
                objectives.append(object)
            pareto_front_solution.append(objectives)
            print(objectives)
            ob=np.asarray(objectives)
            func1_best.append(min(ob[:,0]))
            func2_best.append(min(ob[:,1]))
            # func1_best.append(sum(ob[:,0])//len(ob))
            # func2_best.append(sum(ob[:,1])//len(ob))
            print(pareto_front_solution)
            print(func1_best)
            print(func2_best)
        return pareto_front_solution,func1_best,func2_best


if __name__ == "__main__":
    # files = [ 11, 13, 15]
    # for file_index in range(15,16):
    #     print(f'../1_Brandimarte/BrandimarteMk{file_index}.fjs')
    #     P_time, Job_number, Machine_number, AGV_Trans, AGV_number, AGV_weight, epk, epik, ea_load, ea_unload, ea_idle = Get_Problem(
    #         f'../1_Brandimarte/BrandimarteMk{file_index}.fjs')


    GEN = 400
    ps1 = 100
    pc1 = 0.8
    pm1 = 0.2

    Machine_number = 10
    Job_number = 7
    AGV_number = 3
    AGV_weight = 2
    # epk = [random.randint(10, 18) for i in range(Machine_number)]
    epk = [10, 8, 12, 10, 7, 13, 11, 8, 10, 11]
    # epik = [random.randint(1, 4) for i in range(Machine_number)]
    epik = [2, 1, 2, 1, 2, 1, 1, 2, 3, 2]
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
        [9999, 9999, 4, 6, 9999, 9999, 9999, 9999, 9999, 9999],
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
        [9999, 9999, 7, 6, 7, 9999, 9999, 9999, 9999, 9999],
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





    EH = EHA(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk, epik, ea_load,
            ea_unload, ea_idle)
    pop1, pop2, pop3 = EH.initial_population()  # 初始化
    for p in range(len(pop1)):
        pop1[p] = EH.jbchrom_reset(pop1[p], pop3[p])  # 修整工序染色体
    pt0=[]
    for i in range(len(pop1)):
        pt1 = []
        pt1.append(pop1[i])
        pt1.append(pop2[i])
        pt1.append(pop3[i])
        pt0.append(pt1)
    print(pt0)
    pareto_front_solution,func1_best,func2_best=EH.main(pt0)
    print(pareto_front_solution)
    print(func1_best)
    save_plot_main_data("pic_data/EHA/data/", pareto_front_solution,func1_best,func2_best,GEN)
    #这边代码我自己删除了重复的染色体，导致没有出现需要删除的精英个体，聚类删除部分代码不确定是否正确

