% Function to return, and optionally plot, the 1931 CIE color matching functions.
% Reference: https://en.wikipedia.org/wiki/CIE_1931_color_space#Color_matching_functions
% You can use these curves to build X, Y, and Z images from a collection of individual spectral images.
% With that XYZ image, you can then get the RGB image with the xyz2rgb() function of the Image Processing Toolbox.
% Values taken from:
% http://www.cie.co.at/publ/abst/datatables15_2004/x2.txt
% http://www.cie.co.at/publ/abst/datatables15_2004/y2.txt
% http://www.cie.co.at/publ/abst/datatables15_2004/z2.txt
% With x, y, and z, the first column is wavelength in nm, and the second column is the chromaticity value.
% xyz is a 2D array where each row = nm, x, y, z for each wavelength.
function [x, y, z, nm, xyz] = get_xyz()
try
	% nm, xbar
	x = [...
		380, 0.001368
		385, 0.002236
		390, 0.004243
		395, 0.007650
		400, 0.014310
		405, 0.023190
		410, 0.043510
		415, 0.077630
		420, 0.134380
		425, 0.214770
		430, 0.283900
		435, 0.328500
		440, 0.348280
		445, 0.348060
		450, 0.336200
		455, 0.318700
		460, 0.290800
		465, 0.251100
		470, 0.195360
		475, 0.142100
		480, 0.095640
		485, 0.057950
		490, 0.032010
		495, 0.014700
		500, 0.004900
		505, 0.002400
		510, 0.009300
		515, 0.029100
		520, 0.063270
		525, 0.109600
		530, 0.165500
		535, 0.225750
		540, 0.290400
		545, 0.359700
		550, 0.433450
		555, 0.512050
		560, 0.594500
		565, 0.678400
		570, 0.762100
		575, 0.842500
		580, 0.916300
		585, 0.978600
		590, 1.026300
		595, 1.056700
		600, 1.062200
		605, 1.045600
		610, 1.002600
		615, 0.938400
		620, 0.854450
		625, 0.751400
		630, 0.642400
		635, 0.541900
		640, 0.447900
		645, 0.360800
		650, 0.283500
		655, 0.218700
		660, 0.164900
		665, 0.121200
		670, 0.087400
		675, 0.063600
		680, 0.046770
		685, 0.032900
		690, 0.022700
		695, 0.015840
		700, 0.011359
		705, 0.008111
		710, 0.005790
		715, 0.004109
		720, 0.002899
		725, 0.002049
		730, 0.001440
		735, 0.001000
		740, 0.000690
		745, 0.000476
		750, 0.000332
		755, 0.000235
		760, 0.000166
		765, 0.000117
		770, 0.000083
		775, 0.000059
		780, 0.000042];
	
	% First column is nm, second column is the chromaticity value.
	% nm, ybar
	y = [...
		380, 0.000039
		385, 0.000064
		390, 0.000120
		395, 0.000217
		400, 0.000396
		405, 0.000640
		410, 0.001210
		415, 0.002180
		420, 0.004000
		425, 0.007300
		430, 0.011600
		435, 0.016840
		440, 0.023000
		445, 0.029800
		450, 0.038000
		455, 0.048000
		460, 0.060000
		465, 0.073900
		470, 0.090980
		475, 0.112600
		480, 0.139020
		485, 0.169300
		490, 0.208020
		495, 0.258600
		500, 0.323000
		505, 0.407300
		510, 0.503000
		515, 0.608200
		520, 0.710000
		525, 0.793200
		530, 0.862000
		535, 0.914850
		540, 0.954000
		545, 0.980300
		550, 0.994950
		555, 1.000000
		560, 0.995000
		565, 0.978600
		570, 0.952000
		575, 0.915400
		580, 0.870000
		585, 0.816300
		590, 0.757000
		595, 0.694900
		600, 0.631000
		605, 0.566800
		610, 0.503000
		615, 0.441200
		620, 0.381000
		625, 0.321000
		630, 0.265000
		635, 0.217000
		640, 0.175000
		645, 0.138200
		650, 0.107000
		655, 0.081600
		660, 0.061000
		665, 0.044580
		670, 0.032000
		675, 0.023200
		680, 0.017000
		685, 0.011920
		690, 0.008210
		695, 0.005723
		700, 0.004102
		705, 0.002929
		710, 0.002091
		715, 0.001484
		720, 0.001047
		725, 0.000740
		730, 0.000520
		735, 0.000361
		740, 0.000249
		745, 0.000172
		750, 0.000120
		755, 0.000085
		760, 0.000060
		765, 0.000042
		770, 0.000030
		775, 0.000021
		780, 0.000015];
	
	% First column is nm, second column is the chromaticity value.
	% nm, zbar
	z = [...
		380, 0.006450
		385, 0.010550
		390, 0.020050
		395, 0.036210
		400, 0.067850
		405, 0.110200
		410, 0.207400
		415, 0.371300
		420, 0.645600
		425, 1.039050
		430, 1.385600
		435, 1.622960
		440, 1.747060
		445, 1.782600
		450, 1.772110
		455, 1.744100
		460, 1.669200
		465, 1.528100
		470, 1.287640
		475, 1.041900
		480, 0.812950
		485, 0.616200
		490, 0.465180
		495, 0.353300
		500, 0.272000
		505, 0.212300
		510, 0.158200
		515, 0.111700
		520, 0.078250
		525, 0.057250
		530, 0.042160
		535, 0.029840
		540, 0.020300
		545, 0.013400
		550, 0.008750
		555, 0.005750
		560, 0.003900
		565, 0.002750
		570, 0.002100
		575, 0.001800
		580, 0.001650
		585, 0.001400
		590, 0.001100
		595, 0.001000
		600, 0.000800
		605, 0.000600
		610, 0.000340
		615, 0.000240
		620, 0.000190
		625, 0.000100
		630, 0.000050
		635, 0.000030
		640, 0.000020
		645, 0.000010
		650, 0.000000
		655, 0.000000
		660, 0.000000
		665, 0.000000
		670, 0.000000
		675, 0.000000
		680, 0.000000
		685, 0.000000
		690, 0.000000
		695, 0.000000
		700, 0.000000
		705, 0.000000
		710, 0.000000
		715, 0.000000
		720, 0.000000
		725, 0.000000
		730, 0.000000
		735, 0.000000
		740, 0.000000
		745, 0.000000
		750, 0.000000
		755, 0.000000
		760, 0.000000
		765, 0.000000
		770, 0.000000
		775, 0.000000
		780, 0.000000];
	clc;
	nm = x(:, 1);
	xyz = [x, y(:, 2), z(:, 2)];
	
	% If you want to see the curves plotted, set plotCurves = true.
	% If you DO NOT want to see the curves plotted, set plotCurves = false.
	plotCurves = false;
	if ~plotCurves
		% Exit now if you don't want to see the curves plotted.
		return;
	end
	% Only get here if you want to see the curves plotted.
	fontSize = 25;
	lineWidth = 4;
	plot(nm, x(:, 2), 'r-', 'LineWidth', lineWidth);
	grid on;
	hold on;
	plot(nm, y(:, 2), 'g-', 'LineWidth', lineWidth);
	plot(nm, z(:, 2), 'b-', 'LineWidth', lineWidth);
	title('The 1931 CIE Color Matching Curves', 'FontSize', fontSize);
	xlabel('Wavelength in nanometers', 'FontSize', fontSize);
	ylabel('Chromatic Response', 'FontSize', fontSize);
	legend({'x', 'y', 'z'}, 'FontSize', fontSize, 'FontWeight','bold', 'Location', 'north');
	text(445, 1.7, 'z', 'FontSize', fontSize, 'FontWeight','bold', 'Color', 'b');
	text(550, 1.05, 'y', 'FontSize', fontSize, 'FontWeight','bold', 'Color', 'g');
	text(597, 1.10, 'x', 'FontSize', fontSize, 'FontWeight','bold', 'Color', 'r');
	text(440, 0.4, 'x', 'FontSize', fontSize, 'FontWeight','bold', 'Color', 'r');
	
	% Set up figure properties:
	% Enlarge figure to full screen.
	set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
	% Get rid of tool bar and pulldown menus that are along top of figure.
	set(gcf, 'Toolbar', 'none', 'Menu', 'none');
	% Give a name to the title bar.
	set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')
	
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	uiwait(warndlg(errorMessage));
end

%====================================================================================================
% APPENDIX: Sample code to turn a bunch of separate wavelength images into a true color RGB image.
% 		for j = 1 : totalNumberOfImagesInListBox    % Loop though all selected indexes.
% 			index = Selected(j);    % Get the next selected index.
% 			% Get the filename for this selected index.
% 			baseImageFileName = ListOfImageNames{index};
% 			imageFullFileName = fullfile(handles.imageFolder, baseImageFileName);
% 
% 			% Display the image.
% 			imgOriginal = DisplayImage(handles, imageFullFileName);
% 			% If imgOriginal is empty (couldn't be read), skip to next file.
% 			if isempty(imgOriginal) 
% 				continue;
% 			end
% 
% 			axes(handles.axesImage);
% 			imshow(imgOriginal);
% 			title(baseImageFileName, 'FontSize', fontSize);
% 			
% 			% Find out what wavelength this image is.
% 			this_nm = ParseFileName(imageFullFileName);
% 			wavelengthRow = find(this_nm == nm);
% 			% Find out the x, y, and z weights to use
% 			this_x = x(wavelengthRow, 2);
% 			this_y = y(wavelengthRow, 2);
% 			this_z = z(wavelengthRow, 2);
% 			
% 			% Analyze this image.
% 			if j == 1
% 				XImage = this_x * double(imgOriginal);
% 				YImage = this_y * double(imgOriginal);
% 				ZImage = this_z * double(imgOriginal);
% 			else
% 				XImage = XImage + this_x * double(imgOriginal);
% 				YImage = YImage + this_y * double(imgOriginal);
% 				ZImage = ZImage + this_z * double(imgOriginal);
% 			end
% 			colorXYZImage = cat(3, XImage, YImage, ZImage);
% 			% Normalize the image 
% 			XYZImage = mat2gray(colorXYZImage);
% 			rgbImage = xyz2rgb(XYZImage);
% 			imshow(rgbImage, []);
% 			title('Full Color, RGB Image', 'FontSize', fontSize);
% 			drawnow;
% 		end	