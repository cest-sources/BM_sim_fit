function PLOT_SPACE(Sim,Space,NUMERIC_SPACE,ANALYTIC_SPACE)
% plot function for asym- and z-spectra simulated using BATCH_simulation.m
% last change: 2014/04/03 by PS

assert( isstruct( Sim ), 'Simulation parameter Sim not set.');
assert( ~isempty( fieldnames(NUMERIC_SPACE) ) || ~isempty( fieldnames(ANALYTIC_SPACE) ), 'No data to plot given.' );

% adjust param_xscale for all_offsets = 0 plots
param_xscale={'lin','lin','lin','lin','lin','lin','lin','lin','lin'};


% set FontSize and LineWidth for legend/label/title...
myLineWidth         = 1;
myMarkerSize        = 6;
fontsize_legend     = 12;
fontsize_label      = 14;
fontsize_ticks      = 14;
fontsize_title      = 16;

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXX GET VARIABLES AND LABLES & CREATE COLORMAP  XXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% get variables from Sim struct
all_offsets = Sim.all_offsets;
modelfield  = Sim.modelfield;

% get labels and units for the different parameters
[PLabel, PUnit]=getSimLabelUnit(); 

% get the names and number of elements in Space struct
names = fieldnames(Space);
n_names = numel(names);
n_elements = zeros(1,n_names);

% create color map
for jj=1:n_names
    n_elements(jj) = numel(Space.(names{jj}));  % needed for colormap
end
L = max(n_elements);
cc = lines(L);

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXXXXXX DETERMINE SUBPLOT SIZES  XXXXXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

plot_dim = [round(sqrt(n_names)) ceil(n_names/round(sqrt(n_names)))];

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXX GET SPECTRA FOR ALL PARAMETERS  XXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% parameter loop (e.g. B0, B1, tp...)
for jj=1:n_names
    
    % delete old elements from workspace
    clear Onres zspec ana_z Zsim Zmod h hname hvec hvec2
    % get the current parameter
    field = names{jj}; 
    % get the number of current parameter's values
    n_ii = numel(Space.(field));
    
    % parameter's value loop (e.g. B1 = 0.4, 0.8, 1.6, 3.2...)
    for ii=1:n_ii 
 
        Sim = [];
        xmod = [];
        R1rho = [];
        if isfield( ANALYTIC_SPACE, field )
            Sim     = ANALYTIC_SPACE.(field){ii}.Sim;

            % get analytical Z-spectrum
            xmod  = ANALYTIC_SPACE.(field){ii}.x;
            Zmod  = ANALYTIC_SPACE.(field){ii}.(modelfield);
            try
                R1rho = -ANALYTIC_SPACE.(field){ii}.R1rho;
            catch
                'no R1rho exists'
                R1rho = 1;
            end
        end
        
        xsim = [];
        Zsim = [];
        if isfield( NUMERIC_SPACE, field )
            Sim     = NUMERIC_SPACE.(field){ii}.Sim;
            
            % get numerical Z-spectrum
            xsim  = NUMERIC_SPACE.(field){ii}.x;
            Zsim  = NUMERIC_SPACE.(field){ii}.zspec;
        end
        
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXXXXX PLOT ASYM- AND Z-SPECTRA  XXXXXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        if all_offsets == 1 % if complete z-spectrum
            
            % % % % % % % % %
            % % z-spectra % %
            % % % % % % % % %
            
            % create figure and subplot
            figure(222),
            subplot(plot_dim(1), plot_dim(2), jj),

            % get subplot label for the legend
            if iscell(Sim.(field))
                str          = sprintf('%s = %s %s', PLabel.(field), Sim.(field){1}, PUnit.(field));
            else
                 str          = sprintf('%s = %.2f %s', PLabel.(field), Sim.(field), PUnit.(field));
            
            end;
            hname(ii)    = {str};
            
            % plot numerical AND analytical solution
            if Sim.numeric && Sim.analytic
                % numerical solution
%                 plot(xsim,Zsim,'Marker','d','Markersize', 5,'MarkerFaceColor','black','color',cc(ii,:),'LineStyle','none');
                plot(xsim,Zsim,'Marker','.','Markersize', myMarkerSize,'MarkerFaceColor',cc(ii,:),'color',cc(ii,:),'LineStyle','none');
                hold on; 
                % analytical solution
                h = plot(xmod,Zmod,'LineWidth',myLineWidth,'color',cc(ii,:),'DisplayName',str);
                hvec(ii)= h; % handle vector for the legend
%                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} numerical (\diamondsuit) and analytical (---) solutions'},'FontSize',12);
                
            % plot numerical solution only
            elseif Sim.numeric && ~Sim.analytic
%                 h = plot(xsim,abs(Zsim),'Marker','d','Markersize', 5,'MarkerFaceColor','black','color',cc(ii,:),'LineStyle','-','DisplayName',str);
                h = plot(xsim,abs(Zsim),'Marker','.','Markersize', myMarkerSize,'MarkerFaceColor',cc(ii,:),'color',cc(ii,:),'LineStyle','-','LineWidth',myLineWidth,'DisplayName',str);
                hold on;
                hvec(ii)= h; % handle vector for the legend
%                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} numerical (\diamondsuit) solutions'},'FontSize',12);
                
            % plot analytical solution only
            elseif ~Sim.numeric && Sim.analytic
                h = plot(xmod,Zmod,'LineWidth',myLineWidth,'color',cc(ii,:),'DisplayName',str);
                hold on;
                hvec(ii)= h; % handle vector for the legend
%                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} analytical (---) solutions'},'FontSize',12);
            end
            % general plot settings
            xlabel('\Delta\omega [ppm]','FontSize',fontsize_label);
            ylabel('Z(\Delta\omega)','FontSize',fontsize_label);
            ylim([0 1.1]);
            set(gca,'xdir','reverse');
            set(gca, 'FontSize', fontsize_ticks);

            

            % % % % % % % % %
            % % MTR-asym  % %
            % % % % % % % % %
            
            % create figure and subplot
            figure(223),
            subplot(plot_dim(1), plot_dim(2), jj),
            
            % calculate asym-spectra
            [ax, az]=asym_PS(xsim,Zsim);
            [anax, anaz]=asym_PS(xmod,Zmod);
            
            % plot numerical AND analytical solution
            if Sim.numeric && Sim.analytic
                % numerical solution
%                 plot(ax,az,'Marker','d','Markersize', 5,'MarkerFaceColor','black','color',cc(ii,:),'LineStyle','none');
                plot(ax,az,'.-','Markersize', myMarkerSize,'MarkerFaceColor',cc(ii,:),'color',cc(ii,:),'LineStyle','none');
                hold on;
                % analytical solution
                h2 = plot(anax,anaz,'LineWidth',myLineWidth,'color',cc(ii,:),'DisplayName',str);
                hvec2(ii)= h2; % handle vector for the legend
%                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} numerical (\diamondsuit) and analytical (---) solutions'},'FontSize',12);
             
            % plot numerical solution only  
            elseif Sim.numeric && ~Sim.analytic
%                 h2 = plot(ax,az,'Marker','d','Markersize', 5,'MarkerFaceColor','black','color',cc(ii,:),'LineStyle','none','DisplayName',str);
                h2 = plot(ax,az,'.-','Markersize', myMarkerSize,'MarkerFaceColor',cc(ii,:),'color',cc(ii,:),'LineStyle','-','LineWidth',myLineWidth,'DisplayName',str);
                hold on;
                hvec2(ii)= h2; % handle vector for the legend
%                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} numerical (\diamondsuit) solutions'},'FontSize',12);
                
            % plot analytical solution only
            elseif ~Sim.numeric && Sim.analytic
                h2 = plot(anax,anaz,'LineWidth',myLineWidth,'color',cc(ii,:),'DisplayName',str);
                hold on;
                hvec2(ii)= h2; % handle vector for the legend
%                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} analytical (---) solutions'},'FontSize',12);
            end
            
            % general plot settings
            xlabel('\Delta\omega [ppm]','FontSize',fontsize_label);
            ylabel('MTR_{asym}(\Delta\omega)','FontSize',fontsize_label);
            set(gca,'xdir','reverse');
            set(gca, 'FontSize', fontsize_ticks);
            
            
            
            
            % % % % % % % % %
            % % % AREX  % % %
            % % % % % % % % %
            
            
%             % create figure and subplot
%             figure(224),
%             subplot(plot_dim(1), plot_dim(2), jj),
%             
%             % calculate AREX-spectra
%             [x_arex, y_arex]=AREX_PS(xsim,Zsim,Sim.R1A);
%             [xx_arex, yy_arex]=AREX_PS(xmod,Zmod,Sim.R1A);
%             
%             % plot numerical AND analytical solution
%             if Sim.numeric && Sim.analytic
%                 % numerical solution
% %                 plot(x_arex,y_arex,'Marker','d','Markersize', 5,'MarkerFaceColor','black','color',cc(ii,:),'LineStyle','none');
%                 plot(x_arex,y_arex,'Marker','0','Markersize', 5,'color',cc(ii,:),'LineStyle','none');
%                 hold on;
%                 % analytical solution
%                 h2 = plot(xx_arex,yy_arex,'LineWidth',1.5,'color',cc(ii,:),'DisplayName',str);
%                 hvec2(ii)= h2; % handle vector for the legend
% %                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} numerical (\diamondsuit) and analytical (---) solutions'},'FontSize',12);
%              
%             % plot numerical solution only  
%             elseif Sim.numeric && ~Sim.analytic
% %                 h2 = plot(x_arex,y_arex,'Marker','d','Markersize', 5,'MarkerFaceColor','black','color',cc(ii,:),'LineStyle','none','DisplayName',str);
%                 h2 = plot(x_arex,y_arex,'Marker','.','Markersize', 5,'color',cc(ii,:),'LineStyle','-','LineWidth',1.5,'DisplayName',str);
%                 hold on;
%                 hvec2(ii)= h2; % handle vector for the legend
% %                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} numerical (\diamondsuit) solutions'},'FontSize',12);
%                 
%             % plot analytical solution only
%             elseif ~Sim.numeric && Sim.analytic
%                 h2 = plot(xx_arex,yy_arex,'LineWidth',1.5,'color',cc(ii,:),'DisplayName',str);
%                 hold on;
%                 hvec2(ii)= h2; % handle vector for the legend
% %                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} analytical (---) solutions'},'FontSize',12);
%             end
%             
%             % general plot settings
%             xlabel('\Delta\omega [ppm]','FontSize',11);
%             ylabel('AREX(\Delta\omega)','FontSize',11);
%             set(gca,'xdir','reverse');
            
            % % % % % % % % %
            % % % R1rho % % %
            % % % % % % % % %
            
            if Sim.analytic
                % create figure and subplot
                figure(225),
                subplot(plot_dim(1), plot_dim(2), jj),

                % plot analytical solution only
                h2 = plot(xsim,R1rho,'LineWidth',myLineWidth,'color',cc(ii,:),'DisplayName',str);
                hold on;
                hvec3(ii)= h2; % handle vector for the legend
%                 title({sprintf('%s' ,PLabel.(field));'\fontsize{8} analytical (---) solutions'},'FontSize',12);
                
                % general plot settings
                xlabel('\Delta\omega [ppm]','FontSize',11);
                ylabel('R_{1\rho}(\Delta\omega)','FontSize',11);
                set(gca,'xdir','reverse');
            end
            
            
            
            
        
        elseif all_offsets == 0 % if only one offset
            
            Onres.(field).x(ii)=Sim.(field);
            
            if Sim.numeric
                % Z-spectrum/value
                Onres.(field).Zss(ii,:) = Zsim;
                % ASYM
                Onres.(field).ASYM(ii) = Zsim(2)-Zsim(1);
                % MTRRex
                Onres.(field).MTRRex(ii) = 1/Zsim(1)-1/Zsim(2);
                % AREX
                Onres.(field).AREX(ii) = (1/Zsim(1)-1/Zsim(2)).* Sim.R1A;
                % dAREX
                Onres.(field).dAREX(ii) = Onres.(field).AREX(ii)./(1-exp(-Sim.R1A./(1/(Zsim(1)*Sim.tp*Sim.n)/Sim.DC)));
            end
            
            if Sim.analytic
                % Z-spectrum/value
                Onres.(field).Zss_ana(ii,:) = Zmod;
                % ASYM
                Onres.(field).ASYM_ana(ii) = Zmod(2)-Zmod(1);
                % MTRRex
                Onres.(field).MTRRex_ana(ii) = 1/Zmod(1)-1/Zmod(2);
                % AREX
                Onres.(field).AREX_ana(ii) = (1/Zmod(1)-1/Zmod(2)).* Sim.R1A;
                % dAREX
                Onres.(field).dAREX_ana(ii) = Onres.(field).AREX_ana(ii)./(1-exp(-Sim.R1A./(1/(Zmod(1)*Sim.tp*Sim.n)/Sim.DC))); 
            end
            
        elseif all_offsets == 2 % if on-resonant only (e.g. on-res. Spin-Lock)
            % Z-spectrum/values
            Onres.(field).x(ii)=Sim.(field); 
            Onres.(field).Zss(ii,:)     = Zsim;

            if Sim.numeric
                Onres.(field).Zss(ii,:) = Zsim;
            end
            
            if Sim.analytic
                Onres.(field).Zss_ana(ii,:) = Zmod;
                Onres.(field).R1rho(ii,:) = R1rho; 
            end
            
        end;  % end all_offset cases
        
    end; % end parameter's value loop (ii)
    
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXX ADD LEGENDS XXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    if all_offsets == 1 && plot_dim(1) <= 2 % add legend
        legend(hvec,hname,'Location','SouthEast','FontSize',fontsize_legend);   % legend z-specta
        legend(hvec2,hname,'Location','NorthWest','FontSize',fontsize_legend);  % legend asym-spectra
    end
    
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXX all_offsets = 0 / 2 XXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    if all_offsets == 0

        % plot Z(offset) over varyvals
        figure(333),
        subplot(plot_dim(1), plot_dim(2), jj),
        
        
        if Sim.numeric
        h   = plot(Onres.(field).x,Onres.(field).Zss(:,1),'.','color','blue','DisplayName','numerical solution');
        hold on
        
        
            plot(Onres.(field).x,Onres.(field).Zss(:,2),'.','color','cyan','DisplayName','numerical solution');
        end
        
        if Sim.analytic
            plot(Onres.(field).x,Onres.(field).Zss_ana(:,1),'color','green','DisplayName','analytical solution');
        end
            
        if Sim.numeric && Sim.analytic
            rel1 = (abs((Onres.(field).Zss)))./abs(Onres.(field).Zss_ana)-1;
            ind1 = find( abs(rel1(:,1)) > 0.0001);
            ind2 = find( abs(rel1(:,2)) > 0.0001);
            plot(Onres.(field).x(ind1),abs(rel1(ind1,1)),'rx','LineWidth',myLineWidth,'Markersize',myMarkerSize,'DisplayName','relative deviation');
            plot(Onres.(field).x(ind2),abs(rel1(ind2,2)),'rx','LineWidth',myLineWidth,'Markersize',myMarkerSize);
            
            legend('show','Location','SouthEast');
        end

        xlabel(PLabel.(field));
        ylabel(['Z ( \Delta\omega=' sprintf('%.2f )',Sim.offset) ' ppm']);
        set(gca,'XScale',param_xscale{jj});
        set(gca,'YScale','lin');
        str = sprintf('Z (%s)', PLabel.(field));
        title(str);

        % plot MTR_asym(offset) over varyvals    
        figure(334),
        subplot(plot_dim(1), plot_dim(2), jj),
        
        if Sim.numeric
            plot(Onres.(field).x,Onres.(field).ASYM,'.','color','blue','DisplayName','numerical solution');
            hold on
        end
        
        if Sim.analytic
            plot(Onres.(field).x,Onres.(field).ASYM_ana,'color','green','DisplayName','analytical solution');
            hold on
        end

        if Sim.numeric && Sim.analytic
            rel1 = (abs((Onres.(field).ASYM)))./abs(Onres.(field).ASYM_ana)-1;
            ind1 = find( abs(rel1(:)) > 0.001);
            plot(Onres.(field).x(ind1),abs(rel1(ind1)),'rx','LineWidth',myLineWidth,'Markersize',myMarkerSize,'DisplayName','relative deviation');
            
            legend('show');
        end

        xlabel(PLabel.(field));
        ylabel(['MTR_{asym} ( \Delta\omega=' sprintf('%.2f )',Sim.offset) ' ppm']);
        set(gca,'XScale',param_xscale{jj});
        set(gca,'YScale','lin');
        str = sprintf('MTR_{asym} (%s)', PLabel.(field));
        title(str);
            
        % plot MTRRex(offset) over varyvals    
        figure(335),
        subplot(plot_dim(1), plot_dim(2), jj),
        
        if Sim.numeric
            plot(Onres.(field).x,Onres.(field).MTRRex,'.','color','blue','DisplayName','numerical solution');
            hold on
        end
        
        if Sim.analytic
            plot(Onres.(field).x,Onres.(field).MTRRex_ana,'color','green','DisplayName','analytical solution');
            hold on
        end

        if Sim.numeric && Sim.analytic
            rel1 = (abs((Onres.(field).MTRRex)))./abs(Onres.(field).MTRRex_ana)-1;
            ind1 = find( abs(rel1(:)) > 0.001);
            plot(Onres.(field).x(ind1),abs(rel1(ind1)),'rx','LineWidth',myLineWidth,'Markersize',myMarkerSize,'DisplayName','relative deviation');
            legend('show');
        end

        xlabel(PLabel.(field));
        ylabel('MTR_{R{ex}} ( \Delta\omega_B )');
        set(gca,'XScale',param_xscale{jj});
        set(gca,'YScale','lin');
        str = sprintf('MTR_{R{ex}} (%s)', PLabel.(field));
        title(str);
            
        % plot AREX(offset) over varyvals    
        figure(336),
        subplot(plot_dim(1), plot_dim(2), jj),
        
        if Sim.numeric
            plot(Onres.(field).x,Onres.(field).AREX,'.','color','blue','DisplayName','numerical solution');
            hold on
        end
        
        if Sim.analytic
            plot(Onres.(field).x,Onres.(field).AREX_ana,'color','green','DisplayName','analytical solution');
            hold on
        end

        if Sim.numeric && Sim.analytic
            rel1 = (abs((Onres.(field).AREX)))./abs(Onres.(field).AREX_ana)-1;
            ind1 = find( abs(rel1(:)) > 0.001);
            plot(Onres.(field).x(ind1),abs(rel1(ind1)),'rx','LineWidth',myLineWidth,'Markersize',myMarkerSize,'DisplayName','relative deviation');
            legend('show');
        end

        xlabel(PLabel.(field));
        ylabel('AREX ( \Delta\omega_B )');
        set(gca,'XScale',param_xscale{jj});
        set(gca,'YScale','lin');
        str = sprintf('AREX (%s)', PLabel.(field));
        title(str);
        
    elseif all_offsets == 2
        
        % plot on-resonant Z-value over varyvals
        figure(444),
        subplot(plot_dim(1), plot_dim(2), jj),
        h   = plot(Onres.(field).x,abs(Onres.(field).Zss(:,1)),'.','color','blue','DisplayName','numerical solution');
        hold on
        h2  = plot(Onres.(field).x,abs(Onres.(field).Zss_ana(:,1)),'color','green','DisplayName','analytical solution');

%         rel1 = (abs((Onres.(field).Zss)))./abs(Onres.(field).Zss_ana)-1;
%         ind1 = find( abs(rel1(:,1)) > 0.0001);
%         ind2 = find( abs(rel1(:,2)) > 0.0001);
%         plot(Onres.(field).x(ind1),abs(rel1(ind1,1)),'rx','LineWidth',0.5,'Markersize',5,'DisplayName','relative deviation');
%         plot(Onres.(field).x(ind2),abs(rel1(ind2,2)),'rx','LineWidth',0.5,'Markersize',5);

        xlabel(PLabel.(field));
        ylabel('Z ( \Delta\omega_B )');
        set(gca,'XScale',param_xscale{jj});
        set(gca,'YScale','lin');
        legend('show','Location','SouthEast');
        str = sprintf('Z (%s)', PLabel.(field));
        title(str);
        
        if Sim.analytic
            %plot R1rho
            figure(445),
            subplot(plot_dim(1),plot_dim(2),jj),
            h = plot(Onres.(field).x./4.5e-5, Onres.(field).R1rho,'color','red','DisplayName','R1rho'); %PS: why 4.5e-5?!
            hold on

            xlabel(PLabel.(field));
            ylabel('R_{1\rho}');
            set(gca,'XScale',param_xscale{jj}');
            set(gca,'YScale','lin');
            legend('show','Location','SouthEast');
        end
        

    end; % end cases all_offsets == 0 / 2        
end; % end parameter loop (jj)