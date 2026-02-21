function [lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)
	% given a box to average over, define lat-lon idx min-max
	if(strcmp(box, 'AtlantisII_1km'))
		lat_idx_min = 55;
		lat_idx_max = 155;
		lon_idx_min = 120;
		lon_idx_max = 195;
	elseif(strcmp(box, 'AtlantisII_5km'))
                lat_idx_min = 12;
                lat_idx_max = 33;
                lon_idx_min = 24;
                lon_idx_max = 40;
	elseif(strcmp(box, 'NESM_1km')) % covers all three seamounts
                lat_idx_min = 15;
                lat_idx_max = 200;
                lon_idx_min = 60;
                lon_idx_max = 280;
	elseif(strcmp(box, 'NESM_5km')) % covers all three seamounts
                lat_idx_min = 6;
                lat_idx_max = 41;
                lon_idx_min = 13;
                lon_idx_max = 57;
	elseif(strcmp(box, 'domain_1km')) % covers all three seamounts
                lat_idx_min = 5;
                lat_idx_max = 273;
                lon_idx_min = 5;
                lon_idx_max = 295;
	elseif(strcmp(box, 'domain_5km')) % covers all three seamounts
                lat_idx_min = 3;
                lat_idx_max = 54;
                lon_idx_min = 3;
                lon_idx_max = 58;
	end
end
