function [n,materials]=refractive_index(material,wavelength)

wavelength=(wavelength*1e06);

materials={'Fused silica - FS','N-BK7 - NBK7','Poly(methyl methacrylate) - PMMA'};

if strcmp(material,'FS')
    n=sqrt(1+(0.6961663*wavelength^2)/(wavelength^2-0.0684043^2)+(0.4079426*wavelength^2)/(wavelength^2-0.1162414^2)+(0.8974794*wavelength^2)/(wavelength^2-9.896161^2));
elseif strcmp(material,'NBK7')
    n=sqrt(1+(1.03961212*wavelength^2)/(wavelength^2-0.00600069867)+(0.231792344*wavelength^2)/(wavelength^2-0.0200179144)+(1.01046945*wavelength^2)/(wavelength^2-103.560653));
elseif strcmp(material,'PMMA')
    n=sqrt(1+(0.99654*wavelength^2)/(wavelength^2-0.00787)+(0.18964*wavelength^2)/(wavelength^2-0.02191)+(0.00411*wavelength^2)/(wavelength^2-3.85727));
elseif strcmp(material,'TiO2')
    n=sqrt(5.913 + 0.2441/(wavelength^2 -0.0803));
elseif strcmp(material,'TiO2_PQ')
    n=2.042 + 0.157./(wavelength.^2) -0.009./(wavelength.^4);
elseif strcmp(material,'air')
    n=1;
else
    error('That material is not in the refractive index database, please add it before continuing')
end


end