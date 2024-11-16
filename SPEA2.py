
from Data_transfer import Get_Problem
from FJSP_agv import FJSP_agv,Job,Machine,Agv
import random
import copy
import numpy as np
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
    sorted1 = sorted(I, key=lambda x: values1[x])  
    sorted2 = sorted(I, key=lambda x: values2[x]) 
    distance[0] = 4444444444444444
    distance[len(I) - 1] = 4444444444444444
    for k in range(1, len(I) - 1):
        distance[k] = distance[k] + (values1[sorted1[k + 1]] - values1[sorted1[k - 1]]) / (
                max(values1) - min(values1))
    for k in range(1, len(I) - 1):
        distance[k] = distance[k] + (values2[sorted2[k + 1]] - values2[sorted2[k - 1]]) / (
                max(values2) - min(values2))
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
class SPEA2:
    def __init__(self,Pt,Jn,Mn,At,An,Aw,epk,epik,ea_load,ea_unload,ea_idle,pop_size=200,pc=0.9,pm=0.3):
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
        self.archive_size=100

    def initial_population(self):  # 随机初始化规则
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
        new_pt0=[]
        for i in range(self.Pop_size//2):#遍历一半染色体
            t1 = random.randint(0, len(newpt0) - 1)
            t2 = random.randint(0, len(newpt0) - 1)
            while (t1 == t2):
                t2 = random.randint(0, len(newpt0) - 1)
            new_pop1=[[-1 for n in range(len(newpt0[0][m]))]for m in range(len(newpt0[0]))]
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
                        new_pop1[1][j]=newpt0[t1][1][j]
                        new_pop1[2][j]=newpt0[t1][2][j]
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
                    nested_list_3d=copy.copy(self.machine_use)
                    flattened_list_2d = sp.flatten_3d_to_2d(nested_list_3d)
                    rand_num = random.choice(flattened_list_2d[j])
                    new_pop[i][1][j]=rand_num
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
                if a[agv]+1>AGV_weight:
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
            value=sp.fitness(pop[i])
            func_value1.append(value[0])
            func_value2.append(value[1])
        non_dominated_sorted_solution = fast_non_dominated_sort(func_value1, func_value2)  # 非支配排序
        return non_dominated_sorted_solution,func_value1,func_value2

    def dominated_by(self, func_value1, func_value2):
        if (func_value1[0] <= func_value2[0] and func_value1[1] <= func_value2[1]):
            return 0
        elif (func_value2[0] <= func_value1[0] and func_value2[1] <= func_value1[1]):
            return 1
        else:
            return -1


    def calculate_raw_fitness(self, population):
        strength = [0 for i in range(len(population))]
        raw_fitness=[0 for i in range(len(population))]
        dominations = {}
        func_value = []
        for i in range(len(population)):
            func_value.append( sp.fitness(population[i]))
        for i in range(len(population)):
            dominations[i] = []
            func1 = func_value[i]
            for j in range(len(population)):
                if i != j:
                    func2 = func_value[j]
                    x=sp.dominated_by(func1, func2)
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

    def Mating(self, fit):#锦标赛选择
        Matingpool = []
        num_range = list(range(0, len(fit)))
        while(len(Matingpool)<self.Pop_size//2):
            # 随机选择两个个体，进行比较，好的个体留下,避免重复选择两个个体
            select_numbers = random.sample(num_range, 2)
            p = select_numbers[0]
            q = select_numbers[1]
            if fit[p] < fit[q]:
                best = p
            elif fit[p] > fit[q]:
                best = q
            elif fit[p] > fit[q]:
                best = p
            else:
                best = q
            Matingpool.append(best)
        return Matingpool

    def Pareto_crowd(self, non_dominated_sorted_solution, func1, func2, lr):
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
                range(0, len(non_dominated_sorted_solution[i]))] 
            front22 = sorted(non_dominated_sorted_solution2_1, key=lambda x: crowding_distance_values[i][x])
            front = [non_dominated_sorted_solution[i][front22[j]] for j in
                     range(0, len(non_dominated_sorted_solution[i]))]
            front.reverse()
            for value in front:
                new_solution.append(value)
                if (len(new_solution) == lr):  
                    break
            if (len(new_solution) == lr):
                break
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

    def main(self,pt0,At0):
        if At0==[]:
            non_dominated_sorted_solution,value1,value2 = sp.non_dominate(pt0)
            for i1 in non_dominated_sorted_solution[0]:
                At0.append(pt0[i1])
        else:
            all_pop=pt0+At0
            #删除一样的个体
            all_pop=sp.remove_duplicates(all_pop)
            raw_fitness,func_value=sp.calculate_raw_fitness(all_pop)
            fitness=sp.calculate_fitness_and_distances(all_pop,raw_fitness,func_value)
            all_pop,fitness=sp.sort_population(all_pop,fitness)
            At0=[all_pop[i] for i in range(len(all_pop)) if fitness[i]<1]
            fitness0 = [fitness[i] for i in range(len(all_pop)) if fitness[i] < 1]
            if len(At0)<self.archive_size:
                At0=At0+all_pop[len(At0):self.archive_size]
                fitness0=fitness[:self.archive_size]
            else:
                # At0,fitness0=sp.trim_archive(At0,fitness0)
                At0 = all_pop[0:self.archive_size]
                fitness0 = fitness[0:self.archive_size]
            best_pop=sp.Mating(fitness0)
            pt0=[]
            for b in best_pop:
                pt0.append(At0[b])
        pt01=sp.cross(pt0)
        for i in range(len(pt01)):
            pt0_1=sp.jbchrom_reset(pt01[i][0],pt01[i][2])
            pt01[i][0]=pt0_1
        pt02 = sp.mutation(pt0)  # 变异
        for i in range(len(pt02)):
            pt0_1 = sp.jbchrom_reset(pt02[i][0], pt02[i][2])
            pt02[i][0] = pt0_1
        pt=pt0+pt01+pt02
        pt = sp.remove_duplicates(pt)
        non_dominated_sorted_solution, func1, func2 = sp.non_dominate(pt)
        pt_list = sp.Pareto_crowd(non_dominated_sorted_solution, func1, func2, self.Pop_size)
        pt0 = []
        for l in pt_list:
            pt0.append(pt[l])
        func_value1 = []
        func_value2 = []
        for i in range(len(At0)):
            fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,
                           AGV_weight,
                           epk,
                           epik, ea_load, ea_unload, ea_idle)
            self.Jobs, self.Machines, self.AGVs = fja.decode(At0[i][0], At0[i][1], At0[i][2])
            v1, v2 = fja.total_value()
        return pt0,At0, func_value1, func_value2

