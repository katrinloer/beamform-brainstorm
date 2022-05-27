# beamform-brainstorm
New ideas to address array- and/or source-induced anisotropy in seismic noise beamforming

Beamforming methods have become a popular tool to analyse the ambient seismic noise wavefield. They can be used, for example, to assess the presence and orientation of faults and fractures in the subsurface. This requires estimating anisotropic behaviour of surface waves, that is, their velocities as a function of propagation azimuth. However, the geometry of a seismic array affects estimates of velocity and propagation direction of seismic noise wavefields measured with beamforming techniques, which results in apparent anisotropy estimates. 

The code histosyn2022.m computes synthetic wavefields with isotropic velocities for the stations of a given array and analyses them with a simple beamformer: any anisotropy observed in the synthetics must be an artefact introduced by the geometry of the stations or sources. This apparent anisotropy can then be used to correct anisotropy estimates from real data to obtain a more reliable estimate of actual (structual) anisotropy.
