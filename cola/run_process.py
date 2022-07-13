import sys
import os
import getopt
from random import randint

def getPOI():
    return (1,2)

if __name__=='__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "g:r:o:q:", [])
    except getopt.GetoptError as err:
        print('Usage:')
        sys.exit(2)

    graph_file_name = ''
    output_file_name = 'COLA.output'
    query_file_name = 'query_file_name'
    ratio = 1
    #samples = 1
    #cost_l = 0
    #cost_u = 100
    total_n = 0

    for o, a in opts:
        if o == "-g":
            graph_file_name = a
        elif o == "-r":
            ratio = float(a)
        #elif o == "-n":
        #    samples = int(a)
        #elif o == "-l":
        #    cost_l = int(a)
        #elif o == "-u":
        #    cost_u = int(a)
        elif o == "-o":
            output_file_name = a
        elif o == "-q":
            query_file_name = a
        else:
            assert False, "unhandled option"

    size_limit =10
    json_file_name = graph_file_name #+ '.json'
    #print('./cola/a {} {} {} {}'.format(graph_file_name, json_file_name, 5, 'index'))
    #os.system('./cola/a {} {} {} {}'.format(graph_file_name, json_file_name, 5, 'index'))
    #print('Finish visibility graph construction')

    #json_file_name = graph_file_name

    print('../kl/kl {} {} {}'.format(json_file_name, size_limit, ratio))
    os.system('../kl/kl {} {} {}'.format(json_file_name, size_limit, ratio))
    print('Finish graph partition')
    #wait = input('==========================')

    cola_file_name = graph_file_name + '.cola'
    cut_edge_file_name = graph_file_name + '.cut_edge'
    partion_file_name = graph_file_name + '.partition'
    print('./cola -g {} -gb {} {}'.format(cola_file_name, partion_file_name, cut_edge_file_name))
    os.system('./cola -g {} -gb {} {}'.format(cola_file_name, partion_file_name, cut_edge_file_name))
    print('Finish graph indexing')
    #wait = input('==========================')

    peek = open(cola_file_name, 'r')
    total_n = int(peek.readline())
    #print(total_n)

    co_file_name = graph_file_name + '.co'
    nodemap = open(co_file_name, 'r')
    nodemap.readline()

    nodedict = {}
    while(len(nodedict) < total_n):
        line = nodemap.readline().split()
        x = line[2]
        y = line[3]
        xy = (x,y)
        nodedict[xy] = line[1]

    query_file_name_1 = graph_file_name + '.query1'
    query_file_name_2 = graph_file_name + '.query2'

    queryfile = open(query_file_name_1, 'w')
    queryfile2 = open(query_file_name_2, 'w')

    queryguide = open(query_file_name, 'r')
    total_q = int(queryguide.readline())
    feas_cnt = 0
    for i in range(0, total_q):
        line = queryguide.readline().split()
        xy1 = (line[0], line[1])
        xy2 = (line[2], line[3])
        s = nodedict[xy1]
        t = nodedict[xy2]
        c = line[4]
        l = line[5]
        if(line[6]=='1'):
            queryfile.write('{} {} 0 0 {} {}\n'.format(s, t, c, l))
            feas_cnt += 1
        else:
            queryfile2.write('{} {} 0 0 {} {}\n'.format(s, t, c, l))

    queryfile.close()
    queryfile2.close()

    cola_file_name = graph_file_name + '.cola'
    bg_file_name = 'bg.' + cola_file_name
    cut_edge_file_name = graph_file_name + '.cut_edge'
    partion_file_name = graph_file_name + '.partition'
    index_file_name = cola_file_name + '.index'
    order_file_name = cola_file_name + '.order'

    os.system('./cola -g {} -lb {} {} {} -i {} -o {} -f {} -dq {} {}'.format(cola_file_name, bg_file_name, partion_file_name, cut_edge_file_name, index_file_name, order_file_name, output_file_name, query_file_name_1, feas_cnt))
    os.system('./cola -g {} -lb {} {} {} -i {} -o {} -f {} -dq {} {}'.format(cola_file_name, bg_file_name, partion_file_name, cut_edge_file_name, index_file_name, order_file_name, output_file_name, query_file_name_2, total_q - feas_cnt))
    print('Finish graph query')
    #wait = input('==========================')
