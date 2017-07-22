function [  ] = lnox_amf_paper_figs(  )
% Generate and save all the plots for my lightning AMF effect paper
data_path = fullfile(behr_repo_dir, 'Workspaces','LNOx-AMFs');
save_path = fullfile(PaperPaths.lnox_amf_paper, 'Images');
do_cloud_table = false;
map_fig_2panel = false;

do_tables = true;
do_main = true;
do_supp = true;

% Tabulate 5th and 95th percentile of cloud pressure, surface pressure, and
% modis albedo to justify the ranges of the AMF sensitivity test
if do_cloud_table
    fns = {'GLOBETerpres','CloudPressure','MODISAlbedo'};
    dnums = datenum('2012-06-01'):datenum('2012-06-07');
    results = nan(numel(dnums)+1, 2*numel(fns));
    for a=1:numel(dnums)
        Data = load_behr_file(dnums(a));
        for b=1:numel(fns)
            vals = cat_sat_data(Data, fns{b});
            qs = quantile(vals(:), [0.05, 0.95]);
            b2 = 2*(b-1)+1;
            results(a,b2:b2+1) = qs;
        end
    end
    results(end,:) = mean(results(1:end-1,:),1);
    
    rnames = [cellstr(datestr(dnums, 'yyyy-mm-dd')); {'Mean'}];
    vnames = cell(1, size(results,2));
    vnames(1:2:end) = cellfun(@(s) [s,'_5th'], fns, 'UniformOutput', false);
    vnames(2:2:end) = cellfun(@(s) [s,'_95th'], fns, 'UniformOutput', false);
    
    results = array2table(results,'RowNames', rnames, 'VariableNames', vnames)
end

%%%%%%%%%%%%%%%%%%
%%%%% TABLES %%%%%
%%%%%%%%%%%%%%%%%%
if do_tables
        % The table for differences in AMFs among the DC3 profile and WRF profile
        % matched to the DC3 flights
        amf_files = {'DC3-WRF-0-Hybrid-AMF-Sensitivities-ClearSky.mat',...
                     'DC3-WRF-500-Hybrid-AMF-Sensitivities-ClearSky.mat',...
                     'DC3-WRF-665-Hybrid-AMF-Sensitivities-ClearSky.mat',...
                     'DC3-WRF-500-nudge-Hybrid-AMF-Sensitivities-ClearSky.mat',...
                     'DC3-WRF-500-nudge-2xflash-Hybrid-AMF-Sensitivities-ClearSky.mat'};

        file_legend_id = {'0', '500', '665', '500, nudge', '500, nudge, 2x flashrate'};
        fns = {'amfs_dc3', 'amfs_wrf_ft_dc3_bl', 'amfs_wrf_mt_dc3_bl_ut'};
        col_names = {'Profile', 'Avg. AMF', '$\%\Delta$ AMF vs. DC3', 'Avg. AMF SZA $< 40$\textdegree', '$\%\Delta$ AMF(SZA $< 40$\textdegree) vs. DC3'};
        prof_names = {'DC3','Free Trop. Hybrid-%s', 'Mid. Trop. Hybrid-%s'};

        n_profs = (numel(fns)-1)*numel(amf_files)+1;
        amf_mat = nan(n_profs, 4); % table will have four columns: avg AMF, % diff to DC3, avg AMF SZA < 40, % diff SZA < 40 to DC3
        row_names = cell(n_profs,1);
        row_names{1} = 'DC3';
        ai = 0;
        for f=1:numel(amf_files)
            amfs = load(fullfile(data_path, amf_files{f}));
            if f == 1
                a1 = 1;
            else
                a1 = 2; % skip DC3 profiles in all the later files
            end
            for a=a1:numel(fns)
                ai = ai+1;
                fn_a = fns{a};
                row_names{ai} = sprintf(prof_names{a}, file_legend_id{f});
                
                ss = amfs.(fn_a).SZAs < 40;
                amf_mat(ai,1) = nanmean(amfs.(fn_a).AMFs(:));
                amf_mat(ai,3) = nanmean(amfs.(fn_a).AMFs(ss));
                if a > 1
                    amf_mat(ai,[2,4]) = reldiff(amf_mat(ai,[1,3]), amf_mat(1,[1,3]))*100;
                end
                
            end
        end


        amf_tab_caption = 'Results of the AMF sensitivity tests on the hybrid profiles in Fig. \ref{fig:matched-profiles}';
        amf_tab_label = 'tab:matched-prof-amfs';
        make_latex_table(amf_mat, 'file', fullfile(save_path, '..', 'Laughner-LightningEffectsAMFs.tex'),...
            'colnames', col_names, 'rownames', row_names, 'caption', amf_tab_caption, 'm2l', {'%.2f'},...
            'label', amf_tab_label, 'insert', true, 'mark', 'AMFTABLE', 'lines', {'\tophline', '\middlehline', '\bottomhline'});

%         make_latex_table(amf_mat, 'file', 'test2.txt',...
%             'colnames', col_names, 'rownames', row_names, 'caption', amf_tab_caption,...
%             'label', amf_tab_label, 'insert', false, 'mark', 'AMFTABLE', 'lines', {'\tophline', '\middlehline', '\bottomhline'});
end
%%%%%%%%%%%%%%%%%%%
%%%%% FIGURES %%%%%
%%%%%%%%%%%%%%%%%%%

% Figure 1 - the scatter plot of how much each AMF parameter varies the AMF
% for different profiles
if do_main
        fig1sub(1) = misc_amf_lnox_plots('factor-hists', 'AllWRF-MeanBinnedProfile-AMF-Sensitivities-ClearSky.mat');
        fig1sub(2) = misc_amf_lnox_plots('factor-hists', 'AllWRF-MeanBinnedProfile-AMF-Sensitivities-CloudySky.mat');
        fig1 = combine_plots(fig1sub, 'dims', [1 2]);
        close(fig1sub);
        fig1.Position(3) = 2*fig1.Position(3);
        label_subfigs(fig1,'xshift',0.18);
        save_all_the_formats(fig1, fullfile(save_path, 'ParameterEffect'));
        close(fig1);


        % Figure 2 - the clear and cloudy contour plots for percent difference
        % between 0 and 500 mol/flash and 665/500 mol/flash
        fig2sub = gobjects(1,4);
        [fig2sub(1), fig2sub(2)] = misc_amf_lnox_plots('contours', 'AllWRF-MeanBinnedProfile-AMF-Sensitivities-ClearSky.mat', 'clear', '500-0');
        [fig2sub(3), fig2sub(4)] = misc_amf_lnox_plots('contours', 'AllWRF-MeanBinnedProfile-AMF-Sensitivities-CloudySky.mat', 'cloudy', '500-0');
        fig2 = combine_plots(fig2sub);
        close(fig2sub);
        fig2.Position(3) = fig2.Position(3)*2;
        fig2.Position(4) = fig2.Position(4)*2;
        label_subfigs(fig2, 'xshift', 0.18);
        save_all_the_formats(fig2, fullfile(save_path, 'PerdelAMFContours'));
        close(fig2);
        
        fig2bsub = gobjects(1,4);
        [fig2bsub(1), fig2bsub(2)] = misc_amf_lnox_plots('contours', 'AllWRF-MeanBinnedProfile-AMF-Sensitivities-ClearSky.mat', 'clear', '665-500');
        [fig2bsub(3), fig2bsub(4)] = misc_amf_lnox_plots('contours', 'AllWRF-MeanBinnedProfile-AMF-Sensitivities-CloudySky.mat', 'cloudy', '665-500');
        fig2b = combine_plots(fig2bsub);
        close(fig2bsub);
        fig2b.Position(3) = fig2b.Position(3)*2;
        fig2b.Position(4) = fig2b.Position(4)*2;
        label_subfigs(fig2b, 'xshift', 0.18);
        save_all_the_formats(fig2b, fullfile(save_path, 'PerdelAMFContours-665-500'));
        close(fig2b);

        % Figure 3 - comparing the matched WRF/DC3 profiles and the hybrids
        fig3sub(1) = misc_amf_lnox_plots('no2-profiles');
        xlabel('[NO_2] (pptv)'); % fix the subscript in the x label

        Opts = struct('match_file', 'DC3-Comparison-flashrate-1-molflash-500-iccg-2-newprofile-fixedBC-using_wrfgridcorners.mat',...
            'quantity','NO2','lats_filter','all','strat_filter',false,'fresh_filter', false);
        [profs, iswrf] = misc_wrf_chem_comp_plots('make-hybrid',Opts);
        fig3sub(2) = figure;
        l = gobjects(4,1);
        %profs.wrf_prof(~iswrf.wrf_ft_dc3_bl) = nan;
        %profs.dc3_prof(iswrf.wrf_mt_dc3_bl_ut) = nan;
        l(1) = line(profs.wrf_ft_dc3_bl*1e6, profs.wrf_pres, 'color', [0 0.5 0], 'linewidth', 2,'linestyle','none','marker','^','linewidth',2,'markersize',14);
        l(2) = line(profs.wrf_mt_dc3_bl_ut*1e6, profs.wrf_pres, 'color', 'm', 'linewidth', 2,'linestyle','none','marker','o','linewidth',2,'markersize',10);
        l(3) = line(profs.wrf_prof*1e6, profs.wrf_pres, 'color','r','linewidth',2,'linestyle','-','markersize',10);
        l(4) = line(profs.dc3_prof*1e6, profs.dc3_pres, 'color','b','linewidth',2,'linestyle','-','markersize',10);
        set(gca,'ydir','reverse','fontsize',16);
        xlabel('[NO_2] (pptv)'); ylabel('Pressure (hPa)');
        legend(l, {'Free trop. hybrid', 'Mid-trop. hybrid', 'WRF profile', 'DC3 profile'});
        fig3 = combine_plots(fig3sub, 'dims', [2 1]);
        close(fig3sub);
        fig3.Position(4) = 2*fig3.Position(4);
        label_subfigs(fig3, 'xshift', 0.14);
        save_all_the_formats(fig3, fullfile(save_path, 'ProfileComparison'));
        close(fig3);

        % and just the average WRF profiles
        profs = load(fullfile(data_path, 'Profiles', 'WRF-MeanBinned-Profs_17-22_utc.mat'));
        profs = profs.profs;
        fig = figure;
        l = gobjects(3,1);
        l(1) = line(profs.lightning_000mol*1e6, profs.pressures, 'color', 'k', 'linewidth', 2);
        l(2) = line(profs.lightning_500mol*1e6, profs.pressures, 'color', 'b', 'linewidth', 2, 'linestyle', '--');
        l(3) = line(profs.lightning_665mol*1e6, profs.pressures, 'color', 'r', 'linewidth', 2, 'linestyle', '-.');
        legend(l,{'No LNO_x', '500 mol flash^{-1} LNO_x', '665 mol flash^{-1} LNO_x'});
        set(gca,'ydir','reverse','fontsize',16);
        xlabel('[NO_2] (pptv)');
        ylabel('Pressure (hPa)');
        save_all_the_formats(fig, fullfile(save_path, 'WRF-Mean-Profs'));
        close(fig);
        
        fig = figure;
        fns = {'lightning_000mol', 'lightning_500mol', 'lightning_665mol'};
        for a=1:numel(fns)
            vcd = integPr2(profs.(fns{a}), profs.pressures, profs.pressures(1));
            profs.([fns{a},'_sf']) = profs.(fns{a}) / vcd;
        end
        l = gobjects(3,1);
        l(1) = line(profs.lightning_000mol_sf*1e25, profs.pressures, 'color', 'k', 'linewidth', 2);
        l(2) = line(profs.lightning_500mol_sf*1e25, profs.pressures, 'color', 'b', 'linewidth', 2, 'linestyle', '--');
        l(3) = line(profs.lightning_665mol_sf*1e25, profs.pressures, 'color', 'r', 'linewidth', 2, 'linestyle', '-.');
        legend(l,{'No LNO_x', '500 mol flash^{-1} LNO_x', '665 mol flash^{-1} LNO_x'});
        set(gca,'ydir','reverse','fontsize',16);
        xlabel('Shape factor (unitless)');
        ylabel('Pressure (hPa)');
        save_all_the_formats(fig, fullfile(save_path, 'WRF-Mean-ShapeFactors'));
        close(fig);

        % Figures 4 & 5 - scattering weight vectors
        fig4 = misc_amf_lnox_plots('sw-plots', 'AllWRF-MeanBinnedProfile-AMF-Sensitivities-ClearSky.mat');
        save_all_the_formats(fig4, fullfile(save_path, 'ClearSkySWs'));
        close(fig4);
        fig5 = misc_amf_lnox_plots('sw-plots', 'AllWRF-MeanBinnedProfile-AMF-Sensitivities-CloudySky.mat');
        save_all_the_formats(fig5, fullfile(save_path, 'CloudySWs'));
        close(fig5);

        % Figure 6 - BEHR maps (percent change AMF, absolute change VCD, OMI and
        % MODIS cloud filtering)
        % figs(1) is % diff 500-0, figs(2) is %diff 665-500, figs(3) is abs diff
        % 500-0, figs(4) is abs diff 665 - 0
        amf_keep = false(1,4);
        vcd_keep = false(1,4);
        if map_fig_2panel
            amf_keep(1) = true;
            vcd_keep(3) = true;
            plot_dims = [1 2];
            pos_dims = 3;
        else
            amf_keep(1:2) = true;
            vcd_keep(3:4) = true;
            plot_dims = [2 2];
            pos_dims = 3:4;
        end
        amffigs = misc_amf_lnox_plots('plot-behr-avgs', 'BEHRAMFTrop', 'omi');
        close(amffigs(~amf_keep));
        amffigs = amffigs(amf_keep);

        vcdfigs = misc_amf_lnox_plots('plot-behr-avgs', 'BEHRColumnAmountNO2Trop', 'omi');
        close(vcdfigs(~vcd_keep));
        vcdfigs = vcdfigs(vcd_keep);
        allfigs = [amffigs(:); vcdfigs(:)];
        fig6 = combine_plots(allfigs, 'dims', plot_dims);
        fig6.Position(pos_dims) = 2*fig6.Position(pos_dims);
        colormap(blue_red_cmap);
        close(allfigs);

        % for some reason, calling subplot in this case wipes out the copied
        % object. Instead we need to find the child axes and put them in order top
        % left, top right, bottom left, bottom right
        xx = isgraphics(fig6.Children, 'axes');
        fig6ax = fig6.Children(xx);
        pos = cat(1, fig6ax.Position);
        [~,ord_ind] = sortrows(pos,[-2 1]); % sort by vertical position in reverse order, then by horizontal position in forward order

        for a=1:numel(ord_ind)
            % add the subfigure letters. Make the relative and absolute color
            % ranges the same
            if mod(a,2) == 1
                clim = get(fig6ax(ord_ind(a)),'CLim');
            elseif ~map_fig_2panel
                set(fig6ax(ord_ind(a)),'CLim',clim/4);
            end
            label_axis_with_letter(sprintf('(%s)', char(96+a)), 'xshift', 0.18, 'ax', fig6ax(ord_ind(a)));
        end
        save_all_the_formats(fig6, fullfile(save_path,'OMICloud-Filtered-Diffs'));
        close(fig6)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUPPLEMENTAL FIGS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if do_supp
    % Flash count figure - flash counts with and without nudging (as pcolor)
    % and box plot
    figS2sub = misc_wrf_chem_comp_plots('flashcounts', 'Nudged DC3');
    figS2sub(1).Children(end).Title.String = '';
    figS2sub(2).Children(end).Title.String = '';
    figS2 = combine_plots(figS2sub([1,2,4]), 'dims', [3, 1]);
    %figS2.Position(3) = figS2.Position(3)*1.1;
    figS2.Position(4) = figS2.Position(4)*2;
    label_subfigs( figS2, 'xshift', 0.15 );
    save_all_the_formats(figS2, fullfile(save_path, 'Flashcounts'));
    close(figS2sub);
    close(figS2);

    % Temperature and water vapor mean profiles
    Opts = struct('match_file', 'DC3-Comparison-500mol-nonudge-wTTandQ.mat', 'quantity', 'Temperature', 'lats_filter', 'all',...
        'strat_filter', false, 'fresh_filter', false);
    [~, dc3_temp_means, dc3_temp_pres, ~, wrf_temp_means, wrf_temp_pres] = misc_wrf_chem_comp_plots('grand-profile', Opts);
    Opts.quantity = 'Water vapor mixing ratio';
    [~, dc3_q_means, dc3_q_pres, ~, wrf_q_means, wrf_q_pres] = misc_wrf_chem_comp_plots('grand-profile', Opts);

    Opts.match_file = 'DC3-Comparison-500mol-nudge-2xflashrate-wTTandQ.mat';
    Opts.quantity = 'Temperature';
    [~, ~, ~, ~, wrf_nudge_temp_means, ~] = misc_wrf_chem_comp_plots('grand-profile', Opts);
    Opts.quantity = 'Water vapor mixing ratio';
    [~, ~, ~, ~, wrf_nudge_q_means, ~] = misc_wrf_chem_comp_plots('grand-profile', Opts);

    figS3 = figure;
    subplot(2,1,1);
    l = gobjects(3,1);
    l(1) = line(dc3_temp_means, dc3_temp_pres, 'color', 'b', 'linestyle', '-', 'linewidth', 2);
    l(2) = line(wrf_temp_means, wrf_temp_pres, 'color', 'r', 'linestyle', '-', 'linewidth', 2);
    l(3) = line(wrf_nudge_temp_means, wrf_temp_pres, 'color', 'r', 'linestyle', '--', 'linewidth', 2);
    set(gca,'fontsize',16,'ydir','reverse');
    legend(l, {'DC3', 'WRF - not nudged', 'WRF - nudged'});
    xlabel('Temperature (K)');
    ylabel('Pressure (hPa)');

    subplot(2,1,2);
    l2 = gobjects(3,1);
    l2(1) = line(dc3_q_means, dc3_q_pres, 'color', 'b', 'linestyle', '-', 'linewidth', 2);
    l2(2) = line(wrf_q_means, wrf_q_pres, 'color', 'r', 'linestyle', '-', 'linewidth', 2);
    l2(3) = line(wrf_nudge_q_means, wrf_q_pres, 'color', 'r', 'linestyle', '--', 'linewidth', 2);
    set(gca,'fontsize',16,'ydir','reverse');
    legend(l2, {'DC3', 'WRF - not nudged', 'WRF - nudged'});
    xlabel('Water mixing ratio (kg/kg)');
    ylabel('Pressure (hPa)');

    figS3.Position(4) = figS3.Position(4)*2;
    label_subfigs(figS3, 'xshift', 0.15);
    save_all_the_formats(figS3, fullfile(save_path, 'Temperature-and-water'));
    close(figS3);
end

end

function save_all_the_formats(hfig, filename)
eps_name = [filename, '.eps'];
set(hfig, 'paperpositionmode', 'auto');
print(hfig, '-depsc2', '-loose', filename);

fig_name = [filename, '.fig'];
savefig(hfig, fig_name);

saveas(hfig, filename, 'png');
end
