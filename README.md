**MATdrr**
======

## Description

**MATdrr** is a Matlab package for the damped rank reduction (DRR) method and its several variants. The DRR method has a variety of applications in both exploration and earthquake seismology, including but not limited to seismic denoising, seismic reconstruction, seismic diffraction separation, constrained LSRTM, constrained FWI, etc. The initial version of this package was held at https://github.com/chenyk1990/MATdrr, which is no longer maintained.

The Python counterpart of the package can be found at https://github.com/aaspip/pydrr. 

## Reference
    Huang, W., Wang, R., Chen, Y., Li, H., & Gan, S. (2016). Damped multichannel singular spectrum analysis for 3D random noise attenuation. Geophysics, 81(4), V261-V270.

    Chen, Y., Huang, W., Zhang, D., & Chen, W. (2016). An open-source Matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction. Computers & Geosciences, 95, 59-66.
    
    Chen, Y., Zhang, D., Jin, Z., Chen, X., Zu, S., Huang, W., & Gan, S. (2016). Simultaneous denoising and reconstruction of 5-D seismic data via damped rank-reduction method. Geophysical Journal International, 206(3), 1695-1717.
    
    Chen, Y., Huang, W., Yang, L., Oboue, Y.A.S.I., Saad, O.M., and Chen Y.F. 2023, DRR: An open-source multi-platform package for the damped rank-reduction method and its applications in seismology. Computers & Geosciences, in press.
    
BibTeX:

	@article{huang2016dmssa,
	  title={Damped Multichannel Singular Spectrum Analysis for 3{D} Random Noise Attenuation},
	  author={Weilin Huang and Runqiu Wang and  Yangkang Chen and Huijian Li and Shuwei Gan},
	  journal={Geophysics},
	  volume={81},
	  number={4},
	  issue={4},
	  pages={V261-V270},
	  year={2016},
	  publisher={Society of Exploration Geophysicists}
	}

	@article{chen2016drr5d,
	  title={Simultaneous denoising and reconstruction of 5{D} seismic data via damped rank-reduction method},
	  author={Yangkang Chen and Dong Zhang and Zhaoyu Jin and Xiaohong Chen and Shaohuan Zu and Weilin Huang and Shuwei Gan},
	  journal={Geophysical Journal International},
	  volume={206},
	  number={3},
	  issue={3},
	  pages={1695-1717},
	  year={2016}
	}
	
	@article{chen2016drr3d,
	  title={An open-source Matlab code package for improved rank-reduction 3{D} seismic data denoising and reconstruction},
	  author={Yangkang Chen and Dong Zhang and Weilin Huang and Wei Chen},
	  journal={Computers \& Geosciences},
	  volume={95},
	  pages={59-66},
	  year={2016}
	}

	@article{chen2023drr,
	  title={DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology},
	  author={Yangkang Chen and Weilin Huang and Liuqing Yang and Yapo Abol\'{e} Serge Innocent Obou\'{e} and Omar M. Saad and Yunfeng Chen},
	  journal={Computers \& Geosciences},
	  volume={TBD},
	  pages={in press},
	  year={2023}
	}
-----------
## Copyright
    MATdrr developing team, 2013-present
-----------

## License
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)   

-----------

## Install
Using the latest version

    git clone https://github.com/aaspip/MATdrr
    cd MATdrr
    addpath(genpath('./')); #in Matlab command line
    
-----------
## Examples
    The "demo" directory contains all runable scripts to demonstrate different applications of MATdrr. 

-----------
## Dependence Packages
* Matlab 2015 and later versions

-----------
## Naming criteria
    drr3d_xxx.m corresponds to all functions in 3D.
    drr5d_xxx.m corresponds to all functions in 5D.
    odrr_xxx.m corresponds to all ODRR-related functions.
    drr_xxx.m corresponds to all subroutines.
    localxxx.m corresponds to key functions in local windows.
    test_xxx.m corresponds to all DEMO scripts that are runnable.
    
    
-----------
## Development
    The development team welcomes voluntary contributions from any open-source enthusiast. 
    If you want to make contribution to this project, feel free to contact the development team. 

-----------
## Contact
    Regarding any questions, bugs, developments, collaborations, please contact  
    Yangkang Chen
    chenyk2016@gmail.com

-----------
## Gallery
The gallery figures of the MATdrr package can be found at
    https://github.com/chenyk1990/gallery/tree/main/matdrr
Each figure in the gallery directory corresponds to a DEMO script in the "demo" directory with the exactly the same file name.

The following figure shows a 2D localized denoising example using the MATdrr package. Generated by [demos/test_matdrr_drr2d_win.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr2d_win.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr2d_win.png' alt='comp' width=960/>

The following figure shows a 3D localized denoising example using the MATdrr package. Generated by [demos/test_matdrr_drr3d_win.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr3d_win.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr3d_win.png' alt='comp' width=960/>

The following figure shows a 2D diffraction separation example using the MATdrr package. Generated by [demos/test_matdrr_drr2d_diffraction.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr2d_diffraction.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr2d_diffraction.png' alt='comp' width=960/>

The following figure shows a 3D denoising example using the MATdrr package. Generated by [demos/test_matdrr_drr3d.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr3d.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr3d.png' alt='comp' width=960/>

The following figure shows a 3D denoising and reconstruction example using the MATdrr package. Generated by [demos/test_matdrr_drr3drecon.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr3drecon.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr3drecon.png' alt='comp' width=960/>

The following figure shows a 3D diffraction separation example using the MATdrr package. Generated by [demos/test_matdrr_drr3d_diffraction.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr3d_diffraction.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr3d_diffraction.png' alt='comp' width=960/>

The following figure shows a 3D seismic reconstruction example on the USArray data using the MATdrr package. Generated by [demos/test_matdrr_drr3drecon_usarray.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr3drecon_usarray.m) 
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr3drecon_usarray.png' alt='comp' width=960/>

The following figure shows a 5D denoising example using the MATdrr package. Generated by [demos/test_matdrr_drr5d.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr5d.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr5d.png' alt='comp' width=960/>

The following figure shows a 5D denoising and reconstruction example using the MATdrr package. Generated by [demos/test_matdrr_drr5drecon.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr5drecon.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr5drecon.png' alt='comp' width=960/>

The following figure shows a 3D denoising example using the ODRR method in the MATdrr package. Generated by [demos/test_matdrr_odrr3d.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_odrr3d.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_odrr3d.png' alt='comp' width=960/>

The following figure shows a 3D denoising and reconstruction example using the ODRR method in the MATdrr package. Generated by [demos/test_matdrr_odrr3drecon.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_odrr3drecon.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_odrr3drecon.png' alt='comp' width=960/>

The following figure shows a 5D denoising example using the ODRR method in the MATdrr package. Generated by [demos/test_matdrr_odrr5d.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_odrr5d.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_odrr5d.png' alt='comp' width=960/>

The following figure shows a 5D denoising and reconstruction example using the ODRR method in the MATdrr package. Generated by [demos/test_matdrr_odrr5drecon.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_odrr5drecon.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_odrr5drecon.png' alt='comp' width=960/>

The following figure shows a dealiased interpolation (densification) example using the MATdrr package. Generated by [demos/test_matdrr_drr3drecon_dealiase.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr3drecon_dealiase.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr3drecon_dealiase.png' alt='comp' width=960/>

The following figure shows a deblending example using the localized DRR method in the MATdrr package. Generated by [demos/test_matdrr_drr3d_win_deblending.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr3d_win_deblending.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr3d_win_deblending.png' alt='comp' width=960/>

The following figure shows a DAS denoising example using the localized DRR method in the MATdrr package. Generated by [demos/test_matdrr_drr2d_win_dassafod.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr2d_win_dassafod.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr2d_win_dassafod.png' alt='comp' width=960/>

The following figure shows an off-the-grid (OTG) denoising and reconstruction example on USArray data using the MATdrr package. Generated by [demos/test_matdrr_drr3drecon_otg_usarray.m](https://github.com/aaspip/MATdrr/tree/main/demos/test_matdrr_drr3drecon_otg_usarray.m)
<img src='https://github.com/chenyk1990/gallery/blob/main/matdrr/test_matdrr_drr3drecon_otg_usarray.png' alt='comp' width=960/>

