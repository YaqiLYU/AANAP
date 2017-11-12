# AANAP
An implementation of AANAP in CVPR 2015 paper.
by YaqiLYU

https://www.zhihu.com/question/34535199/answer/135169187

_________________________________________________________________________

Lin C C, Pankanti S U, Ramamurthy K N, et al.
Adaptive as-natural-as-possible image stitching [C]
2015 IEEE Conference on Computer Vision and Pattern Recognition (CVPR). 
IEEE, 2015: 1155-1163.

___________________________________________________________________________

HOW TO START

1. Down APAP code from http://cs.adelaide.edu.au/~tjchin/apap/ 
    PLEASE DOWNLOAD Source codes [MDLT code]
  
2. Upzip mdlt.zip, run main.m, Make sure APAP is working correctly

3. Download Dateset1, Dateset2, Dateset3, unzip three dataset to ./mdlt/images/

3. Copy all *.m to ./mdlt/, convert apap to aanap

4. run AANAP.m, set fast_stitch=1 for speed, set fast_stitch=0 for quality!

Just for fun! 
For commercial purposes, please contact the author.
BTW, Please forgive my poor English!

____________________________________________________________________________

这是CVPR 2015的图像拼接论文AANAP的MATLAB实现，基于APAP代码，论文：

Lin C C, Pankanti S U, Ramamurthy K N, et al. Adaptive as-natural-as-possible image stitching [C]//
2015 IEEE Conference on Computer Vision and Pattern Recognition (CVPR). IEEE, 2015: 1155-1163.

仅供学习和参考，商业用途请联系原论文作者！

by YaqiLYU

我在这里：https://www.zhihu.com/question/34535199/answer/135169187

——————————————————————————————————————————————————————————————————————————

HOW TO START:
1. 下载APAP: 在https://cs.adelaide.edu.au/~tjchin/apap/ 下载APAP的MATLAB代码，解压后得到./mdlt.
    我用的是Source codes [MDLT code]，也就是说，仅支持两图拼接，多图拼接请自行修改，非常简单！
    
2. 请先运行main.m，确保APAP可以正常工作。

3. 下载测试图像库：继续在https://cs.adelaide.edu.au/~tjchin/apap/ 下载三个图像数据库，解压到./mdlt/images目录。

4. 拷贝AANAP的所有.m代码到APAP根目录 ./mdlt，注意，这里的所有matlab文件都拷贝过去。

5. 运行AANAP.m就能得到拼接结果了，我提供了两个选项:
    设置fast_stitch=1，快速拼接模式，效果较差。
    设置fast_stitch=0，高质量拼接模式，速度较慢。
    
玩的开心！
版权归源论文作者所有，请不要用于商业用途。
