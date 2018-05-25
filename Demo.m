clc
close all
clear all

%% ------- BRIQUE 1 --------- %%

brique = rgb2gray(imread('brique.png'));

% Filtre de canny 
% 1. on commence par le filtre gaussien
% 2. on calcule le gradient en x et en y
% 3. on calcule la norme 1 du gradient en x et en y (en valeur absolue)
% 4. on calcule la direction du gradient


% Pour le filtre gaussien :
filtreGauss = fspecial('gaussian', [3 3], 0.5);
img_gauss = filter2(filtreGauss, brique);

% Gradient de l'image filtrée :
[img_dx img_dy] = gradient(img_gauss);

% Intensité du gradient avec une carte des intensités :
Intensite_gradient = abs(img_dx) + abs(img_dy);

% Direction du gradient avec une carte des directions :
eps = 10^(-5);
d = atan(img_dy ./ (img_dx+ eps));
d = d * 180 / pi;
[n,p]=size(d);
for i=1:n
    for j=1:p
        if(d(i,j) < 0)
            d(i,j) = d(i,j) + 360;
        end
    end
end

% Paliers de 45°
d = round(d/45) * 45;

seuilH = 15;
seuilB = 5;

img_Canny = zeros(n,p);

for i=2:n-1
    for j=2:p-1
        switch d(i,j)
            case {0, 180, 360}
                voisin1_i = 0;
                voisin2_i = 0;
                voisin1_j = -1;
                voisin2_j = 1;
            
            case {90, 270}
                voisin1_i = -1;
                voisin2_i = 1;
                voisin1_j = 0;
                voisin2_j = 0;
                
            case {45, 225}
                voisin1_i = 1;
                voisin2_i = -1;
                voisin1_j = -1;
                voisin2_j = 1;
                
            case {135, 315}
                voisin1_i = -1;
                voisin2_i = 1;
                voisin1_j = 1;
                voisin2_j = -1;
        end
        
        voisin1 = Intensite_gradient(i+voisin1_i, j+voisin1_j);
        voisin2 = Intensite_gradient(i+voisin2_i, j+voisin2_j);
        
        if(Intensite_gradient(i,j) > seuilH)
            img_Canny(i,j) = 1;
        end
        
        if(Intensite_gradient(i,j) < seuilB)
            img_Canny(i,j) = 0;
        end 
        
        if(Intensite_gradient(i,j) < seuilH && Intensite_gradient(i,j) > seuilB)
            if(voisin1 == 1 || voisin2 == 1)
                img_Canny(i, j) = 1;
            else
                img_Canny(i, j) = 0;
            end
            
        end
    end
end

%Affichage image et image filtree
figure
subplot(1,2,1)
imshow(brique)
title('Image source')
subplot(1,2,2)
imshow(img_Canny);
title('Image filtree')

%% Resultats Brique 1

GLCM_b = graycomatrix(brique)
carac_b=graycoprops(GLCM_b)
d2=d/45;
GCM_tmp1=cooccurrence (d2, 0, 1, 1);
GCM_Brique1=[GCM_tmp1(1:3,1:3) GCM_tmp1(1:3,7:9);
    GCM_tmp1(7:9,1:3) GCM_tmp1(7:9,7:9)]
carac_Brique1=graycoprops(GCM_Brique1)


%% ------- BRIQUE 2 --------- %%


brique2 = rgb2gray(imread('brique2.jpg'));

% Filtre de canny 
% 1. on commence par le filtre gaussien
% 2. on calcule le gradient en x et en y
% 3. on calcule la norme 1 du gradient en x et en y (en valeur absolue)
% 4. on calcule la direction du gradient


% Pour le filtre gaussien :
filtreGauss = fspecial('gaussian', [3 3], 0.5);
img_gauss = filter2(filtreGauss, brique2);

% Gradient de l'image filtrée :
[img_dx img_dy] = gradient(img_gauss);

% Intensité du gradient avec une carte des intensités :
Intensite_gradient = abs(img_dx) + abs(img_dy);

% Direction du gradient avec une carte des directions :
eps = 10^(-5);
d = atan(img_dy ./ (img_dx+ eps));
d = d * 180 / pi;
[n,p]=size(d);
for i=1:n
    for j=1:p
        if(d(i,j) < 0)
            d(i,j) = d(i,j) + 360;
        end
    end
end

% Paliers de 45°
d = round(d/45) * 45;

seuilH = 15;
seuilB = 5;

img_Canny = zeros(n,p);

for i=2:n-1
    for j=2:p-1
        switch d(i,j)
            case {0, 180, 360}
                voisin1_i = 0;
                voisin2_i = 0;
                voisin1_j = -1;
                voisin2_j = 1;
            
            case {90, 270}
                voisin1_i = -1;
                voisin2_i = 1;
                voisin1_j = 0;
                voisin2_j = 0;
                
            case {45, 225}
                voisin1_i = 1;
                voisin2_i = -1;
                voisin1_j = -1;
                voisin2_j = 1;
                
            case {135, 315}
                voisin1_i = -1;
                voisin2_i = 1;
                voisin1_j = 1;
                voisin2_j = -1;
        end
        
        voisin1 = Intensite_gradient(i+voisin1_i, j+voisin1_j);
        voisin2 = Intensite_gradient(i+voisin2_i, j+voisin2_j);
        
        if(Intensite_gradient(i,j) > seuilH)
            img_Canny(i,j) = 1;
        end
        
        if(Intensite_gradient(i,j) < seuilB)
            img_Canny(i,j) = 0;
        end 
        
        if(Intensite_gradient(i,j) < seuilH && Intensite_gradient(i,j) > seuilB)
            if(voisin1 == 1 || voisin2 == 1)
                img_Canny(i, j) = 1;
            else
                img_Canny(i, j) = 0;
            end
            
        end
    end
end

%Affichage image et image filtree
figure
subplot(1,2,1)
imshow(brique2)
title('Image source')
subplot(1,2,2)
imshow(img_Canny);
title('Image filtree')

%% Resultats Brique 2

GLCM_b2 = graycomatrix(brique2)
carac_b2=graycoprops(GLCM_b2)
d3=d/45;
GCM_tmp2=cooccurrence (d3, 0, 1, 1);
GCM_Brique2=[GCM_tmp2(1:3,1:3) GCM_tmp2(1:3,7:9);
    GCM_tmp2(7:9,1:3) GCM_tmp2(7:9,7:9)]
carac_Brique2=graycoprops(GCM_Brique2)



%% ------- GRAINES --------- %%

graine = rgb2gray(imread('graine.png'));

% Filtre de canny 
% 1. on commence par le filtre gaussien
% 2. on calcule le gradient en x et en y
% 3. on calcule la norme 1 du gradient en x et en y (en valeur absolue)
% 4. on calcule la direction du gradient


% Pour le filtre gaussien :
filtreGauss = fspecial('gaussian', [3 3], 0.5);
img_gauss = filter2(filtreGauss, graine);

% Gradient de l'image filtrée :
[img_dx img_dy] = gradient(img_gauss);

% Intensité du gradient avec une carte des intensités :
Intensite_gradient = abs(img_dx) + abs(img_dy);

% Direction du gradient avec une carte des directions :
eps = 10^(-5);
d = atan(img_dy ./ (img_dx+ eps));
d = d * 180 / pi;
[n,p]=size(d);
for i=1:n
    for j=1:p
        if(d(i,j) < 0)
            d(i,j) = d(i,j) + 360;
        end
    end
end

% Paliers de 45°
d = round(d/45) * 45;

seuilH = 15;
seuilB = 5;

img_Canny = zeros(n,p);

for i=2:n-1
    for j=2:p-1
        switch d(i,j)
            case {0, 180, 360}
                voisin1_i = 0;
                voisin2_i = 0;
                voisin1_j = -1;
                voisin2_j = 1;
            
            case {90, 270}
                voisin1_i = -1;
                voisin2_i = 1;
                voisin1_j = 0;
                voisin2_j = 0;
                
            case {45, 225}
                voisin1_i = 1;
                voisin2_i = -1;
                voisin1_j = -1;
                voisin2_j = 1;
                
            case {135, 315}
                voisin1_i = -1;
                voisin2_i = 1;
                voisin1_j = 1;
                voisin2_j = -1;
        end
        
        voisin1 = Intensite_gradient(i+voisin1_i, j+voisin1_j);
        voisin2 = Intensite_gradient(i+voisin2_i, j+voisin2_j);
        
        if(Intensite_gradient(i,j) > seuilH)
            img_Canny(i,j) = 1;
        end
        
        if(Intensite_gradient(i,j) < seuilB)
            img_Canny(i,j) = 0;
        end 
        
        if(Intensite_gradient(i,j) < seuilH && Intensite_gradient(i,j) > seuilB)
            if(voisin1 == 1 || voisin2 == 1)
                img_Canny(i, j) = 1;
            else
                img_Canny(i, j) = 0;
            end
            
        end
    end
end

%Affichage image et image filtree
figure
subplot(1,2,1)
imshow(graine)
title('Image source')
subplot(1,2,2)
imshow(img_Canny);
title('Image filtree')

%% Resultats Graine :

GLCM_g = graycomatrix(graine)
carac_g=graycoprops(GLCM_g)
d4=d/45;
GCM_tmp3=cooccurrence (d4, 0, 1, 1);
GCM_Graine=[GCM_tmp3(1:3,1:3) GCM_tmp3(1:3,7:9);
    GCM_tmp3(7:9,1:3) GCM_tmp3(7:9,7:9)]
carac_Graine=graycoprops(GCM_Graine)
