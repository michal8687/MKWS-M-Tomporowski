function [Nozzle_Input] = InitialiseNozzleShape(Nozzle_Input)

switch Nozzle_Input.shape_case
    case "conical"
        %Calculating the geometric dimentions of the nozzle curves
        convergent_line_length = sqrt((Nozzle_Input.convergent_length-Nozzle_Input.throat_length)^2 + ((Nozzle_Input.entry_diameter - Nozzle_Input.throat_diameter) / 2 - Nozzle_Input.throat_rounding_radius_1).^2 - Nozzle_Input.throat_rounding_radius_1^2);
        convergent_line_angle = acos((Nozzle_Input.convergent_length-Nozzle_Input.throat_length) / sqrt((Nozzle_Input.convergent_length-Nozzle_Input.throat_length)^2 + ((Nozzle_Input.entry_diameter - Nozzle_Input.throat_diameter) / 2 - Nozzle_Input.throat_rounding_radius_1)^2)) .* ...
            (Nozzle_Input.entry_diameter/2 - Nozzle_Input.throat_diameter/2 - Nozzle_Input.throat_rounding_radius_1) ./ abs(Nozzle_Input.entry_diameter/2 - Nozzle_Input.throat_diameter/2 - Nozzle_Input.throat_rounding_radius_1) ...
            + atan(Nozzle_Input.throat_rounding_radius_1/convergent_line_length);

        divergent_line_length = sqrt(Nozzle_Input.divergent_length^2 + ((Nozzle_Input.exit_diameter - Nozzle_Input.throat_diameter) / 2 - Nozzle_Input.throat_rounding_radius_2).^2 - Nozzle_Input.throat_rounding_radius_2^2);
        divergent_line_angle = acos(Nozzle_Input.divergent_length / sqrt(Nozzle_Input.divergent_length^2 + ((Nozzle_Input.exit_diameter - Nozzle_Input.throat_diameter) / 2 - Nozzle_Input.throat_rounding_radius_2)^2)) .* ...
            (Nozzle_Input.exit_diameter/2 - Nozzle_Input.throat_diameter/2 - Nozzle_Input.throat_rounding_radius_2) ./ abs(Nozzle_Input.exit_diameter/2 - Nozzle_Input.throat_diameter/2 - Nozzle_Input.throat_rounding_radius_2) ...
            + atan(Nozzle_Input.throat_rounding_radius_2/divergent_line_length);

        %Setting the nozzle shape as points of x and r (radial axis)
        Nozzle_Input.points = [linspace(0, convergent_line_length .* abs(cos(convergent_line_angle)), 50)', Nozzle_Input.entry_diameter / 2 - linspace(0, convergent_line_length .* sin(convergent_line_angle), 50)'];%Convergent cone

        Nozzle_Input.points = [Nozzle_Input.points; linspace(convergent_line_length .* abs(cos(convergent_line_angle)), Nozzle_Input.convergent_length - Nozzle_Input.throat_length, 50)', ... 
            Nozzle_Input.throat_diameter .* 0.5 + Nozzle_Input.throat_rounding_radius_1 - sqrt(Nozzle_Input.throat_rounding_radius_1.^2 - (linspace(convergent_line_length .* abs(cos(convergent_line_angle)) + Nozzle_Input.throat_length, Nozzle_Input.convergent_length, 50)' - Nozzle_Input.convergent_length).^2)];%Convergent rounding

        if(Nozzle_Input.throat_length > 0)
            Nozzle_Input.points = [Nozzle_Input.points; linspace(Nozzle_Input.convergent_length - Nozzle_Input.throat_length, Nozzle_Input.convergent_length, 50)', linspace(Nozzle_Input.throat_diameter / 2, Nozzle_Input.throat_diameter / 2, 50)'];%Throat section
        end

        Nozzle_Input.points = [Nozzle_Input.points; linspace(Nozzle_Input.convergent_length, Nozzle_Input.convergent_length + Nozzle_Input.divergent_length - divergent_line_length .* abs(cos(divergent_line_angle)), 50)', Nozzle_Input.throat_diameter .* 0.5 + Nozzle_Input.throat_rounding_radius_2 - ... 
            sqrt(Nozzle_Input.throat_rounding_radius_2.^2 - (linspace(Nozzle_Input.convergent_length, Nozzle_Input.convergent_length + Nozzle_Input.divergent_length - divergent_line_length .* abs(cos(divergent_line_angle)), 50)' - Nozzle_Input.convergent_length).^2)];%Divergent rounding

        Nozzle_Input.points = [Nozzle_Input.points; Nozzle_Input.convergent_length + Nozzle_Input.divergent_length + linspace(-divergent_line_length .* abs(cos(divergent_line_angle)), 0, 50)', ...
            Nozzle_Input.exit_diameter / 2 - linspace(divergent_line_length .* sin(divergent_line_angle), 0, 50)'];%Divergent cone

        %Removing unnecesary points to allow interpolation
        if(Nozzle_Input.throat_length > 0)
            Nozzle_Input.points(250,:) = [];
        end
        Nozzle_Input.points(200,:) = [];
        Nozzle_Input.points(150,:) = [];
        Nozzle_Input.points(100,:) = [];
        Nozzle_Input.points(50,:) = [];
    otherwise
        error("Unknown nozzle shape type");
end
end