%TP serie de fourier
%partie 1 - Calcul Analytique


	n=0:19;
	coef = zeros(length(n(end)+1),1);

	% la fonction de calcul des coef de la serie de Fourier
	coef(n+1) = -(-1).^(n-1)./(pi*(n.^2-1)) - 0.9*(-1).^(n-1)*10./(pi*(n.^2-100)); 
	coef(~isfinite(coef))=1;

	%afficher les resultats
	figure('name', 'Transformee de Fourier');
	figure(1);
	bar(abs(coef));
	title([num2str(n(end)+1), 'premiers coefs de Fourier']);908



%partie 2 - calcul des frequences de la fonction.
	%fonction a analyser
	func = @(x) (cos(2*pi*x) + 0.9*cos(2*pi*10*x));

	%definition du domain d'etude
	x = [0:0.025:10];
	y = func(x);

	%afficher les resultat de la fonction
	figure('name', 'Transformee de Fourier');
	figure(2);
	plot(x,y);
	title('Serie de Fourrier');
	xlabel('x');
	ylabel('y');


%partie 3 - utilisation de la methode de NQUIST pour les frequences
	z = fft(y);
	ech = size(z,2); % echantillonage
	freq_ech = 1/0.025;
	freq_x = [0:ech-1] * (freq_ech / ech);
	freq_y = abs(z/ech);

	%affichage des resultats
	figure('name', 'Transformee de Fourier');
	figure(3);
	plot(freq_x, freq_y);
	title('Frequences');
	xlabel('Frequences Transformee de Fourier');
	ylabel('Amplitude normalisee de la Transformee de Fourier');


%partie 4 - reconstruction de la fonction grace a 'ifft' & 'fft'
	%deltas pour reconstruction
	d = [0.025, 0.05, 0.1]; 

	%affichage des resultats
	figure('name', 'Transformee de Fourier');
	figure(4);

	%loop pour 3 differents deltas de reconstruction
	for i=1:3
		x2 = [0:d(i):10];
		freq_ech = length(x2);
		y2 = func(x2);
		z2 = fft(y2) / freq_ech;
		subplot(3,1,i)
		h = plot(x2, (ifft(z2) * freq_ech));
		set(h, 'color', 'red');
		hold on
		plot(x,y);
		title(['reconstruction de la fonction de base grace a ifft() - Interval de ', num2str(d(i))]);
		legend('signal reconstruit', 'signal de depart');
		xlabel('x');
		ylabel('y');
	endfor




%partie 5 - filtrage des frequences - suppression des 2 pics
	
	%recherche et calcul du premier pic
	z_cg = z;
	z_abs = abs(z_cg);
	ind_1 = find(abs(z_abs) == max(abs(z_abs)));
	
	val_1 = z_cg(ind_1(1)-5:ind_1(1)+5);
	val_abs_1 = z_abs(ind_1(1)-5:ind_1(1)+5);
	
	val_2 = z_cg(ind_1(2)-5:ind_1(2)+5);
	val_abs_2 = z_abs(ind_1(2)-5:ind_1(2)+5);

	z_abs(ind_1(1)-5:ind_1(1)+5) = 0;
	z_cg(ind_1(1)-5:ind_1(1)+5) = 0;

	z_abs(ind_1(2)-5:ind_1(2)+5) = 0;
	z_cg(ind_1(2)-5:ind_1(2)+5) = 0;

	freq_y = z_abs / ech;

	%affichage des resultats sans premier pic
	figure('name', 'Filtrage sans pic #1');
	figure(5);
	plot(freq_x,freq_y);
	title('Filtrage sans pic #1');
	xlabel('Frequence de la Transformee de Fourier');
	ylabel('Amplitude normalisee de la Transformee de Fourier');

	figure('name', 'Reconstruction par ifft sans pic #1');
	figure(6);
	plot(x, ifft(z_cg));
	title('Reconstruction par ifft sans pic #1');
	xlabel('x');
	ylabel('y');

	z_cg(ind_1(1)-5:ind_1(1)+5) = val_1;
	z_cg(ind_1(2)-5:ind_1(2)+5) = val_2;


	%recherche et calcul du deuxieme pic
	ind_2 = find(abs(z_abs) == max(abs(z_abs)));
	
	z_abs(ind_1(1)-5:ind_1(1)+5) = val_abs_1;
	z_abs(ind_1(2)-5:ind_1(2)+5) = val_abs_2;

	val_2 = z_abs(ind_2);

	z_abs(ind_2(1)-5:ind_2(1)+5) = 0;
	z_abs(ind_2(2)-5:ind_2(2)+5) = 0;

	z_cg(ind_2(1)-5:ind_2(1)+5) = 0;
	z_cg(ind_2(2)-5:ind_2(2)+5) = 0;

	freq_y = z_abs / ech;

	%affichage des resultats sans deuxieme pic
	figure('name', 'Reconstruction sans pic #2');
	figure(7);
	plot(freq_x, freq_y);
	title('Reconstruction sans pic #2');
	xlabel('Frequence de la Transformee de Fourier');
	ylabel('Amplitude normalisee de la Transformee de Fourier');

	figure('name', 'Transformee de Fourier');
	figure(8);
	plot(x, ifft(z_cg));
	title('Reconstruction par ifft sans pic #2');
	xlabel('x');
	ylabel('y');

%partie 6 - test avec les donnes fournies 'myData.txt'


	vals = importdata('/home/seventh/all_src/2015-2016/2_maths/mydata.txt');
	vals = vals(1:end-1,:);
	%displaying results
	figure('name', 'transformees de fourier - frequencies stored into mydata.txt');
	figure(9);
  plot(vals(:,1), vals(:,2));
	title('Frequencies stored into mydata.txt');
	xlabel('x');
	ylabel('y');
	%calculation of the frequencies' modules using Nynquist methode
	freq_ech = length(x);
	z_vals = fft(vals(:,2));
	ech = size(vals,1);
	%displaying results
	figure('name', 'transformees de fourier - modules des frequences de la serie');
	figure(10);
	

	%Nynquist formula
	x_scale = 1 / 0.025;
	%calculation of the frequencies' modules using Nynquist methode
	freq_x = [0:ech-1] * x_scale / ech;
	freq_y = abs(z_vals / ech);
	%displaying result
	plot(freq_x, freq_y);
	xlabel('Frequences de la transforme de fourier');
	ylabel('Amplitude de la transforme de fourier');

	title('Modules des frequences');
	%filtering

	filter = find (10 > 1/0.025 * vals(:,1));
	
	z_vals_filter = z_vals;
	z_vals_filter(filter(end): end-filter(end)) = 0;
	filter_y = abs(z_vals_filter / ech);
	%displaying results
	figure('name', 'transformees de fourier - modules de frequences apres application du filtre');
	figure(11);
	
	%Nynquist formula
	x_scale = 1 / 0.025;
	%calculation of the frequencies' modules using Nynquist methode
	freq_x = [0:ech-1] * x_scale / ech;
	freq_y = abs(z_vals_filter / ech);
	%displaying result
	plot(freq_x, freq_y);
	xlabel('Frequences de la transforme de fourier');
	ylabel('Amplitude de la transforme de fourier');


	title('Filtrage sans pic #2');
	%rebuilding fourier series
	figure('name', 'Filtrage sans pic #2');
	figure(12);
	h = plot(vals(:,1), ifft(z_vals_filter));
	set(h,'color', 'red');
	title('Modules des frequences');
	xlabel('x');
	ylabel('y');
	hold on
	plot(vals(:,1), vals(:,2));
	legend('signal filtre', 'signal de depart');