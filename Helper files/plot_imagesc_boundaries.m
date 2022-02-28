function plot_imagesc_boundaries(fname, data, boundaries, cax, lim, cmap, labels, backcolor)
%specify cax, lim and cmap as empty [],  if you do not need them


        fig1=figure();
        fig1.Renderer='Painters';
        boundaries=boundaries-0.5;
        h1=imagesc(data); 
        set(h1, 'AlphaData', 1-isnan(data))
        
        if exist('cax', 'var')
            if ~isempty(cax)
                caxis(cax);
            end
        end
        
                
        if exist('backcolor', 'var')
            if ~isempty(backcolor)
                set(gca, 'Color', backcolor)
            end
        end
        
        if exist('labels', 'var')
            if ~isempty(labels)
                set(gca,'Xtick',1:length(labels), 'Xticklabel', labels);
                xtickangle(90)
                set(gca,'Ytick',1:length(labels), 'Yticklabel', labels);
            else
                set(gca,'Xtick',[],'Ytick',[])
            end
        else
            set(gca,'Xtick',[],'Ytick',[])
        end

        axis('square')
        
        if exist('cmap', 'var')
            if ~isempty(cmap)
                colormap(cmap)
            end
        end
        
        hold on;
        for ii=2:length(boundaries)
            xvar=boundaries(ii-1):0.5:boundaries(ii);
            plot(repmat(boundaries(ii),[1 length(xvar)]), xvar, '-k','LineWidth',2)
            plot(xvar,repmat(boundaries(ii),[1 length(xvar)]),  '-k','LineWidth',2)
        end
        for ii=1:length(boundaries)-1
            xvar=boundaries(ii):0.5:boundaries(ii+1);
            plot(xvar,repmat(boundaries(ii),[1 length(xvar)]),  '-k','LineWidth',2)
            plot(repmat(boundaries(ii),[1 length(xvar)]), xvar,  '-k','LineWidth',2)
        end
        
        if exist( 'lim','var')
            if ~isempty(lim)
                h=fig1;
                set(h.Children, 'Xlim', lim);
                set(h.Children, 'Ylim', lim);
            end
        end
        
        colorbar
        [path,name,ext]=fileparts(fname);
        if strcmp(ext, '.jpg')
            saveas(fig1, fname)
        else
            print(fig1, '-dpdf', fname)
        end