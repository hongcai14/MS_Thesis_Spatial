function[results] = Spatial(data,step,s,i,window_radius,window_min)
%==========================================================================
% Needs the following functions:
% 
%   (1) Featherstone, and 
%   (2) error_calc_fun.
%
%   Input:
%    data := earthquake catalogue (filtered for min and max magnitudes),
%    where column
%
%   (1) := month
%   (2) := day
%   (3) := year
%   (4) := hour 
%   (5) := minute
%   (6) := second
%   (7) := magnitude
%   (8) := latitude
%   (9) := longitude
%   (10 := depth (km), may not necessarily be present (needed)
%
%    step := size (degrees) of latitude-by-longitude squares. 
%    s := a row vector of lengthscales (km).
%    i := indices for s 
%    window_radius := radius (km) of each node.
%    window_min := minimum number of events needed to be in each cell's
%    window radius to continue on with calculation
%
%   Output:
%    results := 
%       results(:,:,1) = b
%       results(:,:,2) = error in b 
%       results(:,:,3) = D2
%       results(:,:,4) = error in D2
%   
%==========================================================================

reference = referenceEllipsoid('WGS84');

lati = 37.0:-step:32.0;
%^node latitudes ('y-values')
loni = -122.0:step:-114.0;
%^node longitudes ('x-values')

results = NaN(length(lati),length(loni),6);
%^prepares...
%   (:,:,1) = b (least squares), 
%   (:,:,2) = b error (least squares),
%   (:,:,3) = D2, 
%   (:,:,4) = D2 error 
%   (:,:,5) = b (maximum likelihood), 
%   (:,:,6) = b error (maximum likelihood)

for ii = 1:length(lati)   
    for jj = 1:length(loni)
        
        distances = NaN(size(data,1),1);
        %^prepares a column vector to save distances (km) between each
        %earthquake's epicenter and current node 
        
        lat1 = lati(ii);
        %^current node's latitude
        lon1 = loni(jj);
        %^current node's longitude
      
        for kk = 1:size(data,1)
            distances(kk) = distance(lat1,lon1,data(kk,8),data(kk,9),reference)/1000;
            %^calculates and assigns the distance (km) between the current
            %node and earthquake kk's epicenter
        end
        
        TEMP = [data distances];
        
        TEMP = TEMP(distances <= window_radius,:);
        %^subset of current data where earthquake epicenters satisfy the
        %specified "window_radius" for the current node
        
        if size(TEMP,1) >= window_min
            %^computes if TEMP contains at least a minimum number of
            %earthquakes, "window_min"
            
            TEMP = sortrows(TEMP,size(TEMP,2));
            %^sorts the rows of TEMP in ascending order based on the elements
            %in the last column 
            
            TEMP = TEMP(1:window_min,1:size(TEMP,2));
             %^restricts TEMP to the hypocenters closest to the current node,
               %and removes last column (distances)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%b-value (Least Squares) calculation%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            M = TEMP(:,7);
            %^magnitudes of subset of current data where earthquake
            %epicenters satisfy the specified "window_radius" for the
            %current node 
            
            Fanta = NaN(length(unique(M)),2);
            Fanta(:,2) = unique(M);
            %^magnitudes are arranged from top to bottom in ascending order
            
            for mm = 1:length(unique(M))
                Fanta(mm,1) = sum(M >= Fanta(mm,2));
            end
            
            b_linear = polyfit(Fanta(:,2),log10(Fanta(:,1)),1);
            b_LS = -b_linear(1);
            
            Sunkist = Fanta';
            
            b_LS_error = error_calc_fun(Sunkist(2,:),log10(Fanta(:,1)));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%D2-value calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [x,y,z] = Featherstone(TEMP(:,8:10),1,2,3);
            X = [x y z];
           
            CI = zeros(length(s),1);
            %^needs to be zeros 
            
            for qq = 1:length(s)
                for ww = 1:size(X,1)
                    for ee = 1:size(X,1)
                        
                        term1 = (X(ww,1) - X(ee,1))^2;
                        term2 = (X(ww,2) - X(ee,2))^2;
                        term3 = (X(ww,3) - X(ee,3))^2;
                        
                        answer = heaviside(s(qq) - ...
                            sqrt(term1 + term2 + term3));
                        
                        CI(qq) = CI(qq) + answer;
                    end
                end
            end
            
            CI = (2/(size(X,1)*(size(X,1)-1)))*CI;
            
            sp = s';
            
            D2_x = log10(sp(i));
            D2_y = log10(CI(i));
            
            D2_linear = polyfit(D2_x,D2_y,1);
            D2 = D2_linear(1);
            %^calculates D2
             
            D2_error = error_calc_fun(log10(s(i)),log10(CI(i)));
            %^calculates error in D2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%b-value (Maximum Likelihood) calculation%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            b_ML = (length(M)*log10(exp(1)))/sum(M - min(M));

            num = (M - mean(M)).^2;
            den = length(M)*(length(M) - 1);

            sigma_squared = sum(num./den);

            b_ML_error = 2.30*b_ML*b_ML*sqrt(sigma_squared);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%Recording%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            results(ii,jj,1) = b_LS;
            results(ii,jj,2) = b_LS_error;
            results(ii,jj,3) = D2;
            results(ii,jj,4) = D2_error;
            results(ii,jj,5) = b_ML;
            results(ii,jj,6) = b_ML_error;
            
        end
    end
end
end