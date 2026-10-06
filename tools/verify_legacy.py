#!/usr/bin/env python3
"""Differentially verify against an actual legacy source snapshot, never archived meshes.

The reference's only algorithm instrumentation is a return AFTER an unchanged
optimization step. Both copies also write observations after each checkpoint.
FLOATING supplies the same std::abs overload to BOTH sources; NATIVE leaves the
compiler's legacy overload resolution intact. No thresholds/geometry are patched.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
BASELINE = 'dfcb3db9b219bf8dbdf4e9130547d7c51f0d63d4'
DRIVER = r'''
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef HEXGEN_FLOATING_ABS
using std::abs;
#endif
class hexGen;
static int observationCount = 0, requestedIterations = 1000;
static void Observe(hexGen&);
#define private public
#include "HexGen.cpp"
#undef private
static void SaveMesh(const Mesh& mesh, int width, const std::string& path) {
    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(&mesh.vNum),sizeof(int));
    out.write(reinterpret_cast<const char*>(&mesh.eNum),sizeof(int));
    out.write(reinterpret_cast<const char*>(mesh.v),mesh.vNum*3*sizeof(double));
    for(int e=0;e<mesh.eNum;++e) out.write(reinterpret_cast<const char*>(mesh.e[e]),width*sizeof(int));
}
static void SaveLeaves(hexGen& gen, const std::string& path) {
    std::ofstream out(path,std::ios::binary);
    out.write(reinterpret_cast<const char*>(&gen.leafNum),sizeof(int));
    for(int i=0;i<gen.leafNum;++i) {
        int id=gen.cutArray[i],level=gen.getLevel[id];
        out.write(reinterpret_cast<const char*>(&id),sizeof(int));
        out.write(reinterpret_cast<const char*>(&level),sizeof(int));
    }
}
static void Observe(hexGen& gen) {
    SaveMesh(gen.octreeMesh,8,"checkpoint-"+std::to_string(++observationCount)+".bin");
    std::ofstream out("threshold-"+std::to_string(observationCount)+".bin",std::ios::binary);
    out.write(reinterpret_cast<const char*>(&ELEM_THRES),sizeof(double));
}
static void Mark(const char* stage, std::chrono::steady_clock::time_point& last) {
    const auto now=std::chrono::steady_clock::now();
    std::cout << "TIMING " << stage << ' ' << std::chrono::duration<double>(now-last).count() << std::endl;
    last=now;
}
static void Projection(hexGen& gen) {
#ifdef REFERENCE_BUILD
    gen.ProjectToIsoSurface("projHex.vtk");
#else
    gen.ProjectToIsoSurface("projHex.vtk",requestedIterations);
#endif
    SaveMesh(gen.octreeMesh,8,"projected.bin");
    const int next=std::rand();
    std::ofstream out("next-random.bin",std::ios::binary);
    out.write(reinterpret_cast<const char*>(&next),sizeof(int));
}
static void BalanceFixture(int seed) {
    hexGen gen(6);
    gen.getLevel.resize(levelId[7]);
    for(int l=1;l<6;++l) std::fill(gen.getLevel.begin()+levelId[l],gen.getLevel.begin()+levelId[l+1],l);
    gen.octreeArray.resize(levelId[7],false);
    std::vector<int> leaves;
    for(int id=0;id<levelId[2];++id) gen.octreeArray[id]=true;
    for(int id=levelId[2];id<levelId[3];++id) leaves.push_back(id);
    std::mt19937 rng(seed);
    for(int round=0;round<18;++round) {
        const int cell=leaves[rng()%leaves.size()],level=gen.getLevel[cell];
        if(level>=5) continue;
        int siblings[8]; gen.RefineBrothers(cell,level,siblings);
        std::unordered_set<int> split(siblings,siblings+8);
        std::vector<int> next;
        for(int id:leaves) {
            if(!split.count(id)) next.push_back(id);
            else for(int c=0;c<8;++c) next.push_back(gen.Child(id,level,c));
        }
        leaves.swap(next);
    }
    // Preserve the legacy padded lists, including stale tails in recursive passes.
    gen.cutArray.assign(65536,0); gen.cutArray1.assign(65536,0);
    std::copy(leaves.begin(),leaves.end(),gen.cutArray.begin()); gen.leafNum=leaves.size();
    SaveLeaves(gen,"balance-before.bin");
    gen.StrongBalancedOctree();
    SaveLeaves(gen,"balance-after.bin");
    gen.BOX_LENGTH_RATIO=100.0/gen.voxelSize;
    gen.START_POINT[0]=gen.START_POINT[1]=gen.START_POINT[2]=0;
    gen.OutputOctree("octree.vtk"); SaveMesh(gen.octreeMesh,8,"octree.bin");
    gen.DualFullHexMeshExtraction("dualFullHex.vtk"); SaveMesh(gen.hexMesh,8,"dualFullHex.bin");
}
static void ProjectionFixture() {
    hexGen gen(6);
    gen.triMesh.Initialize(12,8,3); gen.triMesh.vNum=8; gen.triMesh.eNum=12;
    const double vertices[8][3]={{0,0,0},{100,0,0},{100,100,0},{0,100,0},
        {0,0,100},{100,0,100},{100,100,100},{0,100,100}};
    std::memcpy(gen.triMesh.v,vertices,sizeof(vertices));
    for(int f=0;f<6;++f) for(int t=0;t<2;++t) {
        const int local[2][3]={{0,1,2},{0,2,3}};
        for(int j=0;j<3;++j) gen.triMesh.e[2*f+t][j]=fIdC[f][local[t][j]];
    }
    gen.octreeMesh.Initialize(32,64); gen.octreeMesh.vNum=8;
    gen.octreeMesh.eNum=gen.leafNum=1;
    gen.BOX_LENGTH_RATIO=1; gen.START_POINT[0]=gen.START_POINT[1]=gen.START_POINT[2]=0;
    for(int i=0;i<8;++i) {
        gen.octreeMesh.e[0][i]=i;
        for(int a=0;a<3;++a) gen.octreeMesh.v[i][a]=vertices[i][a];
    }
    std::srand(1); Projection(gen);
}
int main(int argc,char** argv) {
    if(argc!=4 && argc!=5) return 2;
    requestedIterations=std::atoi(argv[3]); std::srand(1);
    std::cout << "ABS_HALF " << abs(0.5) << std::endl;
    if(std::string(argv[1])=="balance") { BalanceFixture(std::atoi(argv[2])); return 0; }
    if(std::string(argv[1])=="projection") { ProjectionFixture(); return 0; }
    hexGen gen(VOXEL_SIZE);
    auto last=std::chrono::steady_clock::now();
    if(std::string(argv[1])=="projection-state") {
        if(argc!=5) return 2;
        gen.ReadRawData(argv[2],"modifiedTri.vtk");
        std::ifstream in(argv[4],std::ios::binary);
        int vertices=0,elements=0;
        in.read(reinterpret_cast<char*>(&vertices),sizeof(int));
        in.read(reinterpret_cast<char*>(&elements),sizeof(int));
        if(!in || vertices<=0 || elements<=0) return 3;
        gen.octreeMesh.Initialize(elements*8,vertices*2);
        gen.octreeMesh.vNum=vertices; gen.octreeMesh.eNum=gen.leafNum=elements;
        in.read(reinterpret_cast<char*>(gen.octreeMesh.v),vertices*3*sizeof(double));
        for(int e=0;e<elements;++e) in.read(reinterpret_cast<char*>(gen.octreeMesh.e[e]),8*sizeof(int));
        if(!in) return 4;
        gen.BOX_LENGTH_RATIO=1;
        SaveMesh(gen.octreeMesh,8,"projection-input.bin");
        // A documented stage replay: both versions start from identical exact
        // mesh bytes and reset the random seed to 1 at the projection boundary.
        std::srand(1); Mark("load_projection_state",last);
        Projection(gen); Mark("projection",last); return 0;
    }
    gen.InitializeOctree(argv[2],"modifiedTri.vtk"); Mark("initialize",last);
    SaveMesh(gen.triMesh,3,"surface.bin");
    { std::ofstream out("curvature.bin",std::ios::binary); out.write(reinterpret_cast<const char*>(gen.triMesh.r),gen.triMesh.vNum*sizeof(double)); }
    gen.GetCellValue(); Mark("refinement",last); SaveLeaves(gen,"leaves-before.bin");
    gen.StrongBalancedOctree(); Mark("balance",last); SaveLeaves(gen,"leaves-after.bin");
    gen.octreeArray.clear(); gen.cutArray1.clear();
    gen.OutputOctree("octree.vtk"); Mark("octree_output",last); SaveMesh(gen.octreeMesh,8,"octree.bin");
    gen.DualFullHexMeshExtraction("dualFullHex.vtk"); Mark("dual_extraction",last); SaveMesh(gen.hexMesh,8,"dualFullHex.bin");
    gen.RemoveOutsideElement("dualHex.vtk"); Mark("remove_outside",last); SaveMesh(gen.octreeMesh,8,"dualHex.bin");
    Projection(gen); Mark("projection",last);
}
'''

def instrument(source, reference, release_grid=False):
    if reference and release_grid:
        # The original grid is dead before the tail-recursive call, but its
        # lifetime retains ~24 GiB per recursion level. Free only that dead
        # allocation so the reference can run on this machine. No decisions,
        # leaf arrays or traversal order are changed.
        needle='\tif (unbalancedNode > 0)\n\t\tStrongBalancedOctree();'
        assert source.count(needle)==1
        source=source.replace(needle,'\tpreVecEightCell.clear();\n'+needle)
    start=source.index('void hexGen::ProjectToIsoSurface(')
    end=source.index('inline void hexGen::InitiateElementValence()',start)
    part=source[start:end]
    needle='\t\t\t\tsmallDist = 114514;\n\t\t\t}'
    assert part.count(needle)==1
    addition='\n\t\t\tObserve(*this);'
    part=part.replace(needle,needle+addition)
    if reference:
        part=part.replace('// modify octreeMesh only', '// modify octreeMesh only\n\tint completedIterations = 0;', 1)
        tail='\t\t}\n\t}\n}'
        assert part.count(tail)==1
        part=part.replace(tail,'\t\t}\n\t\tif (++completedIterations >= requestedIterations) {\n\t\t\toctreeMesh.WriteToVtk("iterationMesh.vtk", BOX_LENGTH_RATIO, START_POINT);\n\t\t\treturn;\n\t\t}\n\t}\n}')
    return source[:start]+part+source[end:]

def build(args, version):
    folder=args.output/'build'/version
    folder.mkdir(parents=True,exist_ok=True)
    files=['HexGen.cpp','HexGen.h','Initialization.h','StaticVars.h','Mesh.cpp','Mesh.h']
    hashes={}
    for name in files:
        if version=='old':
            data=subprocess.check_output(['git','show',f'{args.baseline}:HybridOctree_Hex/{name}'],cwd=ROOT)
        else: data=(ROOT/'HybridOctree_Hex'/name).read_bytes()
        hashes[name]=hashlib.sha256(data).hexdigest()
        if name=='HexGen.cpp': data=instrument(data.decode(),version=='old',args.reference_release_grid).encode()
        (folder/name).write_bytes(data)
    (folder/'driver.cpp').write_text(DRIVER)
    command=[args.compiler,'-std=c++11','-O3','-DNDEBUG','driver.cpp','Mesh.cpp','-o','driver']
    if version=='old': command.insert(1,'-DREFERENCE_BUILD')
    if args.abs_mode=='FLOATING': command.insert(1,'-DHEXGEN_FLOATING_ABS')
    with (folder/'build.log').open('w') as log:
        subprocess.run(command,cwd=folder,stdout=log,stderr=subprocess.STDOUT,check=True)
    (folder/'build.json').write_text(json.dumps({'command':command,'source_sha256':hashes},indent=2))
    return folder/'driver'

def limits():
    resource.setrlimit(resource.RLIMIT_AS,(36*1024**3,36*1024**3))

def run(args, executables, case, mode, input_value, state=None):
    results={}
    for version,executable in executables.items():
        work=args.output/'runs'/case/version
        if work.exists(): shutil.rmtree(work)
        work.mkdir(parents=True)
        command=['/usr/bin/time','-f','%e %M','-o',str(work/'time.txt'),str(executable),mode,str(input_value),str(args.iterations)]
        if state is not None: command.append(str(state))
        print(json.dumps({'event':'started','case':case,'version':version}),flush=True)
        with (work/'stdout.log').open('w') as out,(work/'stderr.log').open('w') as err:
            process=subprocess.Popen(command,cwd=work,stdout=out,stderr=err,preexec_fn=limits,start_new_session=True)
            (work/'process.json').write_text(json.dumps({'pid':process.pid,'command':command}))
            try: returncode=process.wait(timeout=args.timeout)
            except BaseException:
                os.killpg(process.pid,signal.SIGTERM)
                process.wait()
                raise
        if returncode: raise RuntimeError(f'{case}/{version} failed ({returncode}); see {work}')
        timings={}
        for line in (work/'stdout.log').read_text().splitlines():
            if line.startswith('TIMING '):
                _,stage,seconds=line.split(); timings[stage]=float(seconds)
        measured=(work/'time.txt').read_text().split()
        results[version]={'seconds':float(measured[0]),'peak_rss_kib':int(measured[1]),'stages_seconds':timings,
            'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(work.iterdir()) if p.suffix in ('.bin','.vtk')}}
        print(json.dumps({'event':'finished','case':case,'version':version,'seconds':results[version]['seconds']}),flush=True)
    a,b=results['old']['hashes'],results['new']['hashes']
    mismatches=[name for name in sorted(a.keys()|b.keys()) if a.get(name)!=b.get(name)]
    result={'case':case,'abs_mode':args.abs_mode,'iterations':args.iterations,'byte_identical':not mismatches,'mismatches':mismatches,'versions':results}
    (args.output/(case+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'event':'compared','case':case,'byte_identical':not mismatches,'mismatches':mismatches}),flush=True)
    if mismatches: raise RuntimeError(f'Different legacy outputs: {case}: {mismatches}')
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',default=BASELINE)
    parser.add_argument('--suite',choices=['fixtures','pipeline','projection'],default='fixtures')
    parser.add_argument('--abs-mode',choices=['NATIVE','FLOATING'],default='NATIVE')
    parser.add_argument('--compiler',default=os.environ.get('CXX','g++'))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--inputs',type=Path,nargs='*',default=[])
    parser.add_argument('--iterations',type=int,default=1000)
    parser.add_argument('--state-dir',type=Path,help='Prior pipeline runs directory; projection replay reads case/old/dualHex.bin')
    parser.add_argument('--timeout',type=int,default=1800)
    parser.add_argument('--reference-release-grid',action='store_true',help='Free the reference dense corner grid before tail recursion (allocation lifetime only)')
    args=parser.parse_args(); args.output=args.output.resolve()
    if args.iterations<1: parser.error('--iterations must be positive')
    if args.suite!='fixtures' and not args.inputs: parser.error('this suite requires --inputs')
    if args.suite=='projection' and not args.state_dir: parser.error('projection requires --state-dir')
    args.output.mkdir(parents=True,exist_ok=True)
    executables={version:build(args,version) for version in ['old','new']}
    results=[]
    if args.suite=='fixtures':
        for seed in [1,7,42]: results.append(run(args,executables,'balance-'+str(seed),'balance',seed))
        results.append(run(args,executables,'projection','projection','unused'))
        args.iterations += 1
        results.append(run(args,executables,'projection-next-step','projection','unused'))
        args.iterations -= 1
        if results[-1]['versions']['new']['hashes']['projected.bin'] == results[-2]['versions']['new']['hashes']['projected.bin']:
            raise RuntimeError('Projection fixture did not advance by one actual optimization step')
    elif args.suite=='pipeline':
        for path in args.inputs: results.append(run(args,executables,path.stem,'pipeline',path.resolve()))
    else:
        for path in args.inputs:
            state=args.state_dir.resolve()/path.stem/'old'/'dualHex.bin'
            results.append(run(args,executables,path.stem,'projection-state',path.resolve(),state))
    report={'baseline':args.baseline,'abs_mode':args.abs_mode,'compiler':subprocess.check_output([args.compiler,'--version'],text=True).splitlines()[0],
            'reference_release_grid':args.reference_release_grid,
            'suite':args.suite,
            'projection_replay': 'Identical pre-projection binary mesh; both versions reset seed to 1 at this stage.' if args.suite=='projection' else None,
            'method':'Same compiler, flags, seed (1), input, constants and numerical overload. Original algorithms; reference adds only observations and a fixed-step return. Output hashes include ordered leaves and unrounded binary mesh state.',
            'results':results}
    (args.output/'results.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__': main()
