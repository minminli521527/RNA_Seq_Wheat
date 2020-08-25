### (1) 提取序列ID名称(去重复)
gene_list = [] #设置空列表用于存储蛋白名称
with open(r'/home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/HANDS/20200627/0_data/Triticum_aestivum.IWGSC.cdna.B_ID','r') as a:  #读取含有蛋白名称的文件
    for i in a: #for循环逐行遍历
        if not i.startswith('#'): #跳过文件中的注释行
            gene_id = i.split(' ')[0] #按空格将每行信息分开并选取第一列蛋白名称
            if gene_id not in gene_list: #如果临时列表中没有当前元素则追加
                gene_list.append(gene_id) #append函数用于将符合条件的蛋白名称添加到列表
print(len(gene_list)) #查看蛋白数目

### (2) 从fa文件提取对应蛋白的序列
final_seq = {} #设置一个空字典用于存放结果
parse_check = False #设置变量用于判断读取的行是否为所需蛋白名称所在行，首先设为flase
with open(r'/home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/HANDS/20200627/0_data/Triticum_aestivum.IWGSC.cdna1.all.fa','r') as a:
    for i in a:
        if i.startswith('>'): #判断是否以>开头
            seq_ID = i.split(' ')[0][1:] #读取fa文件的蛋白名
            if seq_ID in gene_list: #判断fa文件中的蛋白名称是否在所需蛋白列表中
                parse_check = True #如果包含在gene_list中则将parse_check重新赋值为True
                seq=[] #设置空列表用于储存序列
                final_seq[seq_ID] = seq #设置字典中key对应的value为列表形式
            else:
                parse_check = False #如果fa的蛋白不在gene_list中，parse仍为false
        elif parse_check == True: #如果parse为true（即fa蛋白在gene_list中）
            seq_part= i.split('\n')[0] #将序列添加到序列列表seq中
            seq.append(seq_part)

#将字典中的结果输出到results.fa中
outfile = open('Triticum_aestivum.IWGSC.cdna.B.fa', 'w')

for id in final_seq:
    outfile.write('>%s\n' % id) #设置字符串格式提取id信息（即蛋白名称）
    sequence = final_seq[id] #获取id对应的序列
    for i in sequence:
        outfile.write('%s' % (i)) #获取序列信息并设置格式输出
    outfile.write('\n')
outfile.close()
