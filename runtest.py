#! /usr/bin/env python3
import subprocess
import io



def record(cmd, filename):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    outs, errs = proc.communicate()
    proc.wait()
    # print('errors = {0}\n'.format(errs))
    # print('output = {0}\n'.format(outs))

    outfile = open(filename, 'w')

    ostr = outs.decode("utf-8")
    outfile.write(ostr)


    # print(ostr.split('\n'))



    # buf = io.StringIO(str(outs))

    # print(buf.readline())
    # for line in buf.read():
    #     print(line)
    # outfile.write(str(outs))


def test(cmd, filename):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    outs, errs = proc.communicate()
    proc.wait()


    prc_out = outs.decode("utf-8").split('\n')
    # ref_out = infile.read()


    ref_out = open(filename, 'r')
    for lr,lo in zip(ref_out, prc_out):
        ref = [x.strip() for x in lr[:-1].split()]
        res = [x.strip() for x in lo.split()]
        index = 0
        for x,y in zip(res, ref):
            if x != y and index != 7 and index != 3:
                print('Discrepancy:\n{0}\n{1}'.format(lr[:-1],lo))
                exit(1)
            index += 1

    print('CMD {}: ok'.format(cmd))

        # print('{0} / {1}'.format(lr[:-1],lo))

    # infile = open(filename, 'r')
    # for lo,li in zip(str(outs), infile.read()):
    #     print('{0} / {1}'.format(lo,li))



if __name__ == '__main__':
    # record(['build/bibd'], 'reference/bibd.out')
    # record(['build/tempo_scheduler', 'data/tsptw/tsptw_145.txt', '--input-format tsptw'], 'reference/tsptw.out')
    # record(['build/tempo_scheduler', 'data/jsp/la19.txt', '--input-format', 'jsp'], 'reference/jsp.out')
    # record(['build/tempo_scheduler', 'data/jstl/la16_0_1', '--input-format', 'jstl'], 'reference/jstl.out')
    # record(['build/tempo_scheduler', 'data/path/cr_16_1.txt', '--input-format', 'path'], 'reference/path.out')
    # record(['build/tempo_team_orienteering', 'data/tsptwopt/tsptw_2.txt', '--no-transitivity', '--no-edge-finding'], 'reference/tsptwopt.out')
    # record(['build/tempo_sat', 'data/cnf/uf250-01.cnf', '--verbosity', '3'], 'reference/sat.out')
    # record(['build/tempo_rcpsp', 'data/rcpsp/j30/j309_5.sm'], 'reference/rcpsp.out')
    # record(['build/tempo_scheduler', 'data/osp/hurley/j6-per0-0.txt'], 'reference/osp.out')

    test(['build/bibd'], 'reference/bibd.out')
    test(['build/tempo_scheduler', 'data/tsptw/tsptw_145.txt', '--input-format tsptw'], 'reference/tsptw.out')
    test(['build/tempo_scheduler', 'data/jsp/la19.txt', '--input-format', 'jsp'], 'reference/jsp.out')
    test(['build/tempo_scheduler', 'data/jstl/la16_0_1', '--input-format', 'jstl'], 'reference/jstl.out')
    test(['build/tempo_scheduler', 'data/path/cr_16_1.txt', '--input-format', 'path'], 'reference/path.out')
    test(['build/tempo_team_orienteering', 'data/tsptwopt/tsptw_2.txt', '--no-transitivity', '--no-edge-finding'], 'reference/tsptwopt.out')
    test(['build/tempo_sat', 'data/cnf/uf250-01.cnf', '--verbosity', '3'], 'reference/sat.out')
    test(['build/tempo_rcpsp', 'data/rcpsp/j30/j309_5.sm'], 'reference/rcpsp.out')
    test(['build/tempo_scheduler', 'data/osp/hurley/j6-per0-0.txt'], 'reference/osp.out')
    

