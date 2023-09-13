% Copyright (c) 2021 iCoSMoS Lab UMD
% Author : Nakul Garg

function optimizeStencil()
    close all; clear all; clc;
    
    number_of_patterns = 5;
    number_of_holes = 20;

    % Select the metric to optimize, amplitude vs phase
    global type_of_pattern;
    type_of_pattern = "amplitude";
    % type_of_pattern = "phase";

    freqList = [80:20:7920];

    if(1)
        [vals, angles, holes_positions] = get_point_cloud(freqList, number_of_holes, number_of_patterns);
        % disp("Saving patterns...");
        % save(path_to_save);
    else
        % disp("Loading saved patterns...");
        % load(path_to_load)
    end

    disp("Calculating Metrics...");
    M = getMetric(vals, number_of_patterns, angles)./100;
    [max_M, max_M_id] = max(M)

    if(0)   %Metric values for all hole paterns
        figure; hold on;
        FontSize = 20;
        set(gca, 'FontSize', FontSize);
        xlabel("Hole Pattern ID");
        ylabel("Metric Value");
        title("Metric Plot")
        plot(M);
        scatter(max_M_id, max_M);
    end
    if(1)
        for i = [1:length(M)]
            figure; hold on;
            FontSize = 26;
            set(gca, 'FontSize', FontSize);
            xlabel("Frequency (Hz)");
            % ylabel("Phase");
            % ylim([-3, 3]);
            % ylabel("Amplitude");
            % ylim([0, 4]);
            % title(strcat("H_ID = ", num2str(i), ", Metric = ", num2str(M(i))), 'Interpreter', 'none');
            PlotAngles = 120:20:180;
            % PlotAngles = 1:60:length(angles)-170;
            % colorlist = linspecer(length(PlotAngles), 'qua');
            legendnumber = 1;
            for a = PlotAngles
                plot(freqList, squeeze(vals(i, a, :)), 'LineWidth', 2);
                legendstring{legendnumber} = num2str(angles(a));
                legendnumber = legendnumber + 1;
            end
            % legend(legendstring, 'location','best');
            legend(["0 deg", "20 deg", "40 deg", "60 deg"], 'location','best');
            if(1)
                show_hole_pattern(holes_positions(i,:,:));
                title(strcat("H_ID = ", num2str(i), ", Metric = ", num2str(M(i))), 'Interpreter', 'none');
            end
        end
    end
end

function M = getMetric(vals, number_of_patterns, angles)
    M = [];
    for i = 1 :number_of_patterns
        M(i) = 0;
        for a = 1:length(angles)-1
            vals_one = squeeze(vals(i,a+1,:)) - mean(squeeze(vals(i,a+1,:)));
            vals_one = vals_one./std(vals_one);
            vals_two = squeeze(vals(i,a,:)) - mean(squeeze(vals(i,a,:)));
            vals_two = vals_two./std(vals_two);
            M(i) = M(i) + norm(vals_one - vals_two);
            % M(i) = M(i) + norm(squeeze(vals(i,a+1,:)) - squeeze(vals(i,a,:)));
        end
    end
end

function [vals, angles, holes_positions] = get_point_cloud(freqList, number_of_holes, number_of_patterns)
    seed_number = 4321;
    rng(seed_number);

    % cylinder_inner_radius = 0.025;
    % cylinder_inner_radius = 0.02;
    cylinder_inner_radius = 0.0175;
    % cylinder_inner_radius = 0.015;
    cylinder_outer_radius = 0.03;
    cylinder_height = 0.06;
    hole_width = 0.002;

    frequency = freqList;
    lambdas = 340./frequency;

    angles = 0:1:360;

    holes_positions = [];
    vals = [];

    cardoid_amplitudes = get_cardoid();

    % for i = 1:1000
    %     get_hole_pattern(number_of_holes, cylinder_height);
    % end

    for i = 1:number_of_patterns
        disp("Calculating Pattern " + num2str(i));
        if(i==1)
            holes_positions(i,:,:) = get_hole_pattern(number_of_holes, cylinder_height, 11);
        elseif(i==2)
            holes_positions(i,:,:) = get_hole_pattern(number_of_holes, cylinder_height, 356);
        else
            holes_positions(i,:,:) = get_hole_pattern(number_of_holes, cylinder_height, 0);
        end
        vals(i,:,:) = get_frequency_pattern(squeeze(holes_positions(i,:,:)), cylinder_inner_radius, cylinder_outer_radius, lambdas, angles, cardoid_amplitudes);
    end
end

function [vals] = get_frequency_pattern(holes_positions, cylinder_inner_radius, cylinder_outer_radius, lambdas, angles, cardoid_amplitudes)
    vals = [];
    for a = 1:length(angles)
        new_holes_positions = rotate_shell(holes_positions, angles(a));
        for i = 1:length(lambdas)
            vals(a, i) = get_diffraction_pattern(new_holes_positions, cylinder_inner_radius, cylinder_outer_radius, lambdas(i), cardoid_amplitudes);
        end
        % vals(a,:) = vals(a,:) - mean(vals(a,:));
        vals(a,:) = vals(a,:)./std(vals(a,:));
    end
end

function Ifinal = get_diffraction_pattern(holes_positions, cylinder_inner_radius, cylinder_outer_radius, lambda, cardoid_amplitudes)
    D = cylinder_inner_radius;
    wavelength = lambda;

    global type_of_pattern;

    integral1D = 0;
    for i = 1:length(holes_positions)
        hole_pos = holes_positions(i,1);
        hole_angle = holes_positions(i,2);
        r = sqrt((D^2) + (hole_pos^2));
        if(hole_angle<=90 || hole_angle>=270)
            added_distance = cylinder_outer_radius*(1-cos(deg2rad(hole_angle)));
        else
            if(hole_angle<=180)
                added_distance = cylinder_outer_radius + (cylinder_outer_radius*deg2rad(hole_angle-90));
            else
                added_distance = cylinder_outer_radius + (cylinder_outer_radius*deg2rad(abs(270-hole_angle)));
            end
            % added_distance = cylinder_outer_radius*deg2rad(hole_angle);
        end
        hole_contribution = (D/(1j*wavelength*(r^2))) .* exp(-1j*2*pi/wavelength*(r+added_distance));
        hole_contribution = hole_contribution*cardoid_amplitudes(hole_angle+1);
        integral1D = integral1D + hole_contribution;
    end
    if(type_of_pattern == "amplitude")
        Ifinal = abs(integral1D);
    else
        Ifinal = angle(integral1D);
    end
end

function new_holes_positions = rotate_shell(holes_positions, a)
    new_holes_positions = holes_positions;
    new_holes_positions(:,2) = mod(new_holes_positions(:,2)+a, 360);
end

function holes_positions = get_hole_pattern(number_of_holes, cylinder_height, skip_patterns)
    holes_positions = [];
    % for skip_number = 1:skip_patterns
    %     random_angles = floor(360.*rand(1,number_of_holes));
    %     random_heights = (cylinder_height*0.1) + (cylinder_height*0.8).*rand(1,number_of_holes);
    % end
    random_angles = floor(360.*rand(1,number_of_holes));
    random_heights = (cylinder_height*0.1) + (cylinder_height*0.8).*rand(1,number_of_holes);
    holes_positions = [random_heights', random_angles'];

    if(0)
        figure;
        grid = zeros(360,360);
        y_coordinates = floor(360.*holes_positions(:,1)./cylinder_height);
        x_coordinates = holes_positions(:,2);
        for i = 1:number_of_holes
            grid(x_coordinates(i), y_coordinates(i)) = 1;
        end
        imagesc(grid);
    end
end

function show_hole_pattern(holes_positions)
    cylinder_inner_radius = 0.015;
    cylinder_outer_radius = 2*cylinder_inner_radius;
    cylinder_height = 0.06;

    [X,Y,Z] = cylinder(cylinder_outer_radius, 360);
    Z = Z*cylinder_height;

    figure; hold on;
    FontSize = 20;
    set(gca, 'FontSize', FontSize);
    testsubject = surf(X,Y,Z); 
    set(testsubject,'FaceAlpha',0.2);
    set(testsubject,'EdgeColor','none');

    hole_heights = squeeze(holes_positions(:,:,1));
    hole_angles = squeeze(holes_positions(:,:,2));
    scatter3(X(1,hole_angles), Y(1,hole_angles), hole_heights, 'filled', 'MarkerEdgeColor', '#FFFFFF', 'MarkerFaceColor', 'w');
    view(45,10);
end

function cardoid_amplitudes = get_cardoid()
    load('./cardoid_used.mat', 'cardoid_amplitude');
    cardoid_amplitudes = cardoid_amplitude;
end