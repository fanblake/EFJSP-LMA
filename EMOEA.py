#He,论文中出现了OBL策略
"""A multiobjective evolutionary algorithm for achieving energy
efficiency in production environments integrated with multiple
automated guided vehicles
策略主要有两个：
编码采用随机密钥的方式，因此对应的采用了OBL策略
其次采用了SPEA的更新精英的决策方式"""
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
import math

# 拥挤度选择
def crowding_distance(values1, values2, I,distance_func1,distance_func2):
    distance = [0 for i in range(0, len(I))]
    # sorted1 = sort_by_values(I, values1[:])
    # logger.info(f"sorted1:{sorted1}")
    # sorted2 = sort_by_values(I, values2[:])
    # logger.info(f"sorted2:{sorted2}")
    sorted1 = sorted(I, key=lambda x: values1[x])  
    sorted2 = sorted(I, key=lambda x: values2[x])   
    # logger.info(f"sorted1_new:{sorted1_new}")
    # logger.info(f"sorted1_new:{sorted2_new}")
    distance[0] = 4444444444444444
    distance[len(I) - 1] = 4444444444444444
    for k in range(1, len(I) - 1):
        distance[k] = distance[k] + (values1[sorted1[k + 1]] - values1[sorted1[k - 1]]) / (
                distance_func1)
    for k in range(1, len(I) - 1):
        distance[k] = distance[k] + (values2[sorted2[k + 1]] - values2[sorted2[k - 1]]) / (
                distance_func2)#这边我进行了修整，不知道对不对，跟参考的代码编写不同，我感觉我的是对的
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
from Data_transfer import Get_Problem
from FJSP_agv import FJSP_agv,Job,Machine,Agv
import random
import copy
import numpy as np
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

class EMOEA:
    def __init__(self,Pt,Jn,Mn,At,An,Aw,epk,epik,ea_load,ea_unload,ea_idle,pop_size=200,pc=0.9,pm=0.2):
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

    def initial_population(self):  # 随机初始化规则
        Pop1 = []
        Pop2 = []
        Pop3 = []
        # 工件的初始化规则：
        for i in range(self.Pop_size):
            random_num = [random.random() for _ in range(sum(self.op_num))]  # 生成随机数
            sorted_indices = sorted(range(len(random_num)), key=lambda i: random_num[i])
            sorted_indices = [sorted_indices[i] + 1 for i in range(sum(self.op_num))]
            new_M = []
            max_num = max(self.op_num)
            for i in range(1, max_num + 1):
                for j, num in enumerate(self.op_num):
                    if i <= num:
                        new_M.append(j)
            new_M.pop(0)
            new_M.append(0)  # 生成这样的形式[2, 3, 4, 1, 2, 3, 4, 1, 3, 1]-1
            chrom_0 = [-1 for i in range(len(sorted_indices))]
            s = 0
            while (s < len(sorted_indices)):
                index = sorted_indices.index(s+1)
                chrom_0[index] = new_M[s]
                s += 1
            chrom = [x for x in chrom_0 for _ in range(2)]  # 染色体生成
            Pop1.append(chrom)
            machine_order = []
            for j in range(len(self.machine_use)):
                for k in range(len(self.machine_use[j])):
                    rand_num = random.choice(self.machine_use[j][k])
                    machine_order.append(rand_num)
            Pop2.append(machine_order)
            AGV_order = [random.randint(0, self.An - 1) for i in range(sum(self.op_num))]
            Pop3.append(AGV_order)
        return Pop1, Pop2, Pop3

    def cross(self, At0):
        new_pt0 = []
        for i in range(self.Pop_size // 2):  # 遍历一半染色体
            t1 = random.randint(0, len(At0) - 1)
            t2 = random.randint(0, len(At0) - 1)
            while (t1 == t2):
                t2 = random.randint(0, len(At0) - 1)
            new_pop1 = [[-1 for n in range(len(At0[0][m]))] for m in range(len(At0[0]))]
            new_pop2 = [[-1 for n in range(len(At0[0][m]))] for m in range(len(At0[0]))]
            gj = [i for i in range(self.fja.Jn)]
            if np.random.rand() < self.pc:
                random.shuffle(gj)  # MPX交叉
                # 工序交叉
                split_point = random.randint(1, len(gj) - 1)
                x1 = gj[:split_point]
                x2 = gj[split_point:]
                for j in range(len(At0[t1][0])):
                    if At0[t1][0][j] in x1:
                        new_pop1[0][j] = At0[t1][0][j]
                site = [j for j, var in enumerate(new_pop1[0]) if var == -1]
                e = 0
                for j in range(len(At0[t1][0])):
                    if At0[t2][0][j] in x2:
                        new_pop1[0][site[e]] = At0[t2][0][j]
                        e = e + 1
                for j in range(len(At0[t1][0])):
                    if pt0[t2][0][j] in x2:
                        new_pop2[0][j] = At0[t2][0][j]
                site = [j for j, var in enumerate(new_pop2[0]) if var == -1]
                e = 0
                for j in range(len(At0[t1][0])):
                    if At0[t1][0][j] in x1:
                        new_pop2[0][site[e]] = At0[t1][0][j]
                        e = e + 1
                # 机器与AGV的交叉
                random_list = [random.choice([0, 1]) for _ in range(sum(self.op_num))]
                for j in range(len(random_list)):
                    if random_list[j] == 1:
                        new_pop1[1][j] = At0[t1][1][j]
                        new_pop1[2][j] = At0[t1][2][j]
                        new_pop2[1][j] = At0[t2][1][j]
                        new_pop2[2][j] = At0[t2][2][j]
                    else:
                        new_pop1[1][j] = At0[t2][1][j]
                        new_pop1[2][j] = At0[t2][2][j]
                        new_pop2[1][j] = At0[t1][1][j]
                        new_pop2[2][j] = At0[t1][2][j]
                new_pt0.append(new_pop1)
                new_pt0.append(new_pop2)
            else:
                new_pt0.append(At0[t1])
                new_pt0.append(At0[t2])
        return new_pt0

    def flatten_3d_to_2d(self, lst):
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
                    flattened_list_2d = EM.flatten_3d_to_2d(nested_list_3d)
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
                            # logger.info("调整染色体")
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
        fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk,
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
            value=EM.fitness(pop[i])
            func_value1.append(value[0])
            func_value2.append(value[1])
        non_dominated_sorted_solution = fast_non_dominated_sort(func_value1, func_value2)  # 非支配排序
        return non_dominated_sorted_solution,func_value1,func_value2

    def dominated_by(self, func_value1, func_value2):
        if (func_value1[0]<=func_value2[0] and func_value1[1]<=func_value2[1]):
            return 0
        elif(func_value2[0]<=func_value1[0] and func_value2[1]<=func_value1[1]):
            return 1
        else:
            return -1

    def calculate_raw_fitness(self, population):
        strength = [0 for i in range(len(population))]
        raw_fitness=[0 for i in range(len(population))]
        dominations = {}
        func_value = []
        for i in range(len(population)):
            func_value.append( EM.fitness(population[i]))

        for i in range(len(population)):
            dominations[i] = []
            func1 =  func_value[i]
            for j in range(len(population)):
                if i != j:
                    func2 = func_value[j]
                    x=EM.dominated_by(func1, func2)
                    if x==0:
                        strength[i] = strength[i] + 1#受i支配的解的数量
                    elif x==1 :
                        dominations[i].append(j)#j支配i
                    else:
                        pass
        for i in range(len(population)):
            for j in dominations[i]:
                raw_fitness[i] = raw_fitness[i] + strength[j]#R(i)等于支配该个体的所有个体的强度值之和
        return raw_fitness,func_value

    def calculate_fitness_and_distances(self, population, raw_fitness,func_value):
        k_n = int(len(population) ** 0.5) - 1
        fitness = [0 for i in range(len(population))]
        all_distances = []
        for i in range(len(population)):
            distances = []
            for j in range(len(population)):
                if i != j:
                    func1 = func_value[i]
                    func2 = func_value[j]

                    distance = sum([(func1[i] - func2[i]) ** 2 for i in
                                    range(2)]) ** 0.5
                    distances.append(distance)
            distances.sort()
            all_distances.append(distances)
            fitness[i] = raw_fitness[i] + 1 / (distances[k_n] + 2)
        return fitness

    def sort_population(self, population, fitness):
        sorted_ids = np.argsort(fitness)
        new_fitness = [fitness[sorted_ids[i]] for i in range(len(population))]
        new_population = [population[sorted_ids[i]] for i in range(len(population))]
        return new_population, new_fitness

    # def trim_archive(self, archive, fitness):#应该是采用聚类方法消除多余的
    #     while len(archive) > self.archive_size:
    #         all_distances = []
    #         for i in range(len(archive)):
    #             distances = []
    #             for j in range(len(archive)):
    #                 if i != j:
    #                     func1 = EM.fitness(archive[i])
    #                     func2 = EM.fitness(archive[j])
    #                     distance = sum([(func1[i] - func2[i]) ** 2 for i in
    #                                     range(2)])
    #                     distances.append(distance)
    #             distances.sort()
    #             all_distances.append(distances)
    #         k_n = 1
    #         while k_n<=len(all_distances)-1:
    #             closest_n = min(all_distances[i][k_n] for i in range(len(all_distances)))
    #             most_crowded = [i for i in range(len(all_distances)) if all_distances[i][k_n] == closest_n]
    #             if len(most_crowded) == 1:
    #                 archive.pop(most_crowded[0])
    #                 fitness.pop(most_crowded[0])
    #                 break
    #             else:
    #                 k_n = k_n + 1
    #     return archive, fitness
    def Mating(self, fit):#锦标赛选择
        Matingpool=[]
        num_range = list(range(0, len(fit)))
        selected_pairs = []  # 用于存储已经被选择的个体
        while(len(num_range)>1):
            #随机选择两个个体，进行比较，好的个体留下,避免重复选择两个个体
            select_numbers=random.sample(num_range,2)
            selected_pairs.append(select_numbers)
            for num in select_numbers:
                num_range.remove(num)
            p=select_numbers[0]
            q=select_numbers[1]
            if fit[p] < fit[q]:
                best = p
            else:
                best = q
            Matingpool.append(best)
        return Matingpool

    def Pareto_crowd(self,non_dominated_sorted_solution,func1,func2,lr):
        crowding_distance_values = []
        logger.info(f"计算距离begin {len(non_dominated_sorted_solution)}")
        distance_func1 = max(func1[:]) - min(func1[:])
        distance_func2 = max(func2[:]) - min(func2[:])
 
        for j in range(0, len(non_dominated_sorted_solution)):
            # logger.info(f"计算距离 param: {func1[:]}, {func2[:]}, {non_dominated_sorted_solution[j][:]}")
            crowding_distance_value = crowding_distance(func1[:], func2[:], non_dominated_sorted_solution[j][:],distance_func1,distance_func2)
            crowding_distance_values.append(crowding_distance_value)  # 拥挤度计算
        logger.info("计算距离end")    
        # 存储层级rank与拥挤度crowd
        rank = [-1 for i in range(len(func1))]
        crowd = [-1 for i in range(len(func1))]
        for j in range(0, len(non_dominated_sorted_solution)):
            for k in range(0, len(non_dominated_sorted_solution[j])):
                rank[non_dominated_sorted_solution[j][k]] = j
                crowd[non_dominated_sorted_solution[j][k]] = crowding_distance_values[j][k]
        new_solution = []
        logger.info("计算距离end1")    
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
                if (len(new_solution) == lr):  # 其实就是根据排序和拥挤度取到足够的个体
                    break
            if (len(new_solution) == lr):
                break
        logger.info("计算距离end2")   
        return new_solution
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
    def main(self,pt0,At0,pareto_front_solution, func1_best, func2_best):
        pt1 = EM.cross(pt0)#交叉
        for i in range(len(pt1)):
            pt0_1 = EM.jbchrom_reset(pt1[i][0], pt1[i][2])
            pt1[i][0] = pt0_1
        pt2 = EM.mutation(pt0)#变异
        for i in range(len(pt2)):
            pt0_1 = EM.jbchrom_reset(pt2[i][0], pt2[i][2])
            pt2[i][0] = pt0_1
        pt=pt1+pt2
        #首先非支配排序选出10%的个体
        non_dominated_sorted_solution,func1,func2= EM.non_dominate(pt)
        #10%的个体进行OBL策略
        logger.info("10%的个体进行OBL策略1")
        new_list=EM.Pareto_crowd(non_dominated_sorted_solution,func1,func2,len(pt)*0.1)
        logger.info("10%的个体进行OBL策略2")
        new_pop=[]
        for l in new_list:
            new_pop.append(pt[l])#需要进行OBL策略的个体
        for p in range(len(new_pop)):
            new_pop1=[]
            new_pop1.append(new_pop[p][0][::-1])#没有采用论文里面的策略，因为很复杂
            new_pop1.append(new_pop[p][1])
            new_pop1.append(new_pop[p][2][::-1])
            pt0_1 = EM.jbchrom_reset(new_pop1[0], new_pop1[2])#修复染色体
            new_pop1[0] = pt0_1
            func_value1 = []
            func_value2 = []
            value = EM.fitness(new_pop[p])
            func_value1.append(value[0])
            func_value2.append(value[1])
            value = EM.fitness(new_pop1)
            func_value1.append(value[0])
            func_value2.append(value[1])
            value=EM.dominated_by(func_value1, func_value2)
            if value==0:
                pass#保持染色体不变
            elif value==1:
                new_pop[p][0]=new_pop1[0]#更换染色体
            else:
                pass#保持染色体不变
        logger.info("10%的个体进行OBL策略 end")
        #进行精英选择策略
        pt=pt+new_pop
        logger.info("去除重复值")
        # pt = EM.remove_duplicates(pt)
        logger.info("non_dominate")
        non_dominated_sorted_solution,func1,func2 = EM.non_dominate(pt)
        logger.info("Pareto_crowd")
        pt_list=EM.Pareto_crowd(non_dominated_sorted_solution,func1,func2,self.Pop_size)
        logger.info("Pareto_crowd end")
        pt0=[]
        for l in pt_list:
            pt0.append(pt[l])
        logger.info("精英选择策略")
        if At0 == []:
            non_dominated_sorted_solution, func1,func2 = EM.non_dominate(pt0)
            for i1 in non_dominated_sorted_solution[0]:
                At0.append(pt0[i1])
        else:
            all_pop=pt0+At0
            #删除一样的个体
            all_pop=EM.remove_duplicates(all_pop)
            logger.info("开始计算原始适应度")
            raw_fitness,raw_func=EM.calculate_raw_fitness(all_pop)
            logger.info("开始计算原始适应距离")
            fitness=EM.calculate_fitness_and_distances(all_pop,raw_fitness,raw_func)
            logger.info("排序")
            all_pop,fitness=EM.sort_population(all_pop,fitness)
            logger.info("排序 end")
            At0=[all_pop[i] for i in range(len(all_pop)) if fitness[i]<1]
            fitness0 = [fitness[i] for i in range(len(all_pop)) if fitness[i] < 1]
            if len(At0)<self.archive_size:
                At0=At0+all_pop[len(At0):self.archive_size]
                fitness0=fitness0+fitness[len(At0):self.archive_size]
            else:
                At0 = all_pop[0:self.archive_size]
                fitness0 = fitness[0:self.archive_size]
                # At0,fitness0=EM.trim_archive(At0,fitness0)
        f_value1 = []
        f_value2 = []
        for i in range(len(At0)):
            fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight,
                           epk,
                           epik, ea_load, ea_unload, ea_idle)
            self.Jobs, self.Machines, self.AGVs = fja.decode(At0[i][0], At0[i][1], At0[i][2])
            v1, v2 = fja.total_value()
            f_value1.append(v1)
            f_value2.append(v2)
        non_dominated_sorted_solution = fast_non_dominated_sort(f_value1, f_value2)  # 非支配排序
        non_dominated_solution = []
        for i in non_dominated_sorted_solution[0]:
            pop = []
            pop.append(At0[i][0])
            pop.append(At0[i][1])
            pop.append(At0[i][2])
            non_dominated_solution.append(pop)
        objectives = []
        for j in non_dominated_sorted_solution[0]:
            object = []
            object.append(f_value1[j])
            object.append(f_value2[j])
            objectives.append(object)
        pareto_front_solution.append(objectives)
        ob = np.asarray(objectives)
        func1_best.append(min(ob[:,0]))
        func2_best.append(min(ob[:,1]))
        logger.info(pareto_front_solution)
        logger.info(func1_best)
        logger.info(func2_best)
        return pt0,At0, pareto_front_solution,func1_best,func2_best


if __name__ == "__main__":
    # files = [ 11, 13, 15]
    # for file_index in range(15,16):
    #     logger.info(f'../1_Brandimarte/BrandimarteMk{file_index}.fjs')
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

    EM = EMOEA(P_time, Job_number, Machine_number, AGV_Trans, AGV_number, AGV_weight, epk, epik, ea_load,
            ea_unload, ea_idle)
    pop1, pop2, pop3 = EM.initial_population()  # 初始化
    for p in range(len(pop1)):
        pop1[p] = EM.jbchrom_reset(pop1[p], pop3[p])  # 修整工序染色体
    pt0=[]
    for i in range(len(pop1)):
        pt1 = []
        pt1.append(pop1[i])
        pt1.append(pop2[i])
        pt1.append(pop3[i])
        pt0.append(pt1)
    logger.info(pt0)
    At0=[]
    pareto_front_solution = []
    func1_best = []
    func2_best = []
    for g in range(GEN):
        logger.info(f'这是算法的第几次迭代*****************************************: {g}')
        pt0,At0,pareto_front_solution, func1_best, func2_best=EM.main(pt0,At0,pareto_front_solution, func1_best, func2_best)
    save_plot_main_data("pic_data/EMOEA/data/", pareto_front_solution, func1_best, func2_best,GEN)
    #这边代码我自己删除了重复的染色体，导致没有出现需要删除的精英个体，聚类删除部分代码不确定是否正确,
    # 为啥最后是两个最优的个体

