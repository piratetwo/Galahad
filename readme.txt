matlab2014a安装包见.iso文件 破解文件见crack文件夹

0.1 安装matlab(工具箱全选)
将$GALAHAD/src/matlab   $CUTEST/src/matlab添加到matlab路径

0.2 安装gfortran-4.7和gcc-4.7并设置其为默认编辑器

0.3 archdefs cutest galahad hslarchive-galahad sif sifdecode这几个文件夹需要在同一目录下
(各个打包文件见.zip文件或使用解压后的文件， hsl-galahad解压后的文件要放入到hslarchive-galahad文件夹中)

注意：(现有目录下已经添加过了，如果是用的压缩包解压后需要)
需要将addfiles文件夹中的 fa14d.f fa14s.f la04d.f la04s.f la15d.f la15s.f  mc29d.f mc29s.f 复制到$GALAHAD/arc/dum文件夹下   la04d.f la04s.f这两个文件原始文件夹下有但需要替换
galahad_lpa_matlab_types.F90  galahad_lpa.F90 galahad_lpa.m复制到$GALAHAD/arc/matlab文件夹下


1.设置环境变量，安装后还需要添加其它的
export GALAHAD=[full path to ./cutest directory]
export ARCHDEFS=[full path to ./archdefs directory]
export SIFDECODE=[full path to ./sifdecode directory]
export CUTEST=[full path to ./cutest directory]
export MATLAB=[full path to matlab directory]($MATLAB/bin下包含mex)
export MYMATLAB=[full path to matlab directory] ($MYMATLAB/bin下包含mex)
export MASTSIF=[full path to ./sif directory]


2. cd $GALAHAD
3. $ARCHDEFS/install_optrove 
4. 依次选择 
y    (1)Everything (b) except AMPL (6) generic 64(不确定平台上的环境) (3)Linux y 
(3)R2013b-2016a  gfortran-4.7 y gcc-4.7 y  
compile SIFDecode y 
compile CUTEst y  D double precision version
compile GALAHAD y

5. 按提示添加环境变量

我这里最后添加的环境变量列表
 export ARCHDEFS="/home/pirate/archdefs"
 export SIFDECODE="/home/pirate/sifdecode"
 export CUTEST="/home/pirate/cutest"
 export GALAHAD="/home/pirate/galahad"
 export PATH="${SIFDECODE}/bin:${PATH}"
 export PATH="${CUTEST}/bin:${PATH}"
 export PATH="${GALAHAD}/bin:${PATH}"
 export MANPATH="${SIFDECODE}/man:${MANPATH}"
 export MANPATH="${CUTEST}/man:${MANPATH}"
 export MANPATH="${GALAHAD}/man:${MANPATH}"

 export MASTSIF="/home/pirate/sif"

 export MYARCH="pc64.lnx.gfo"
 export MYMATLABARCH="pc64.lnx.gfo"

 export MYMATLAB="/usr/local/MATLAB/R2014a"
 export MATLAB="/usr/local/MATLAB/R2014a"

 export MATLABPATH="${CUTEST}/src/matlab:$(MATLABPATH)"
 export MATLABPATH="${GALAHAD}/src/matlab:$(MATLABPATH)"


export MATLABPATH="$MASTSIF:$MATLAB"
export MATLABPATH="${CUTEST}/src/matlab:$MATLAB"
export MATLABPATH="${GALAHAD}/src/matlab:$MATLAB"
