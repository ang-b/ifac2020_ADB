%% output
figure(1)
set(1, 'defaultTextInterpreter', 'latex');
set(1, 'Unit', 'inches', 'Position', [0 0 6 4]);
% we plot just subsystem 1 and attacked
f1_subss = [1 nA];
for i=1:length(f1_subss)
    yT1_rs = reshape(yT(1,f1_subss(i),:), [1 nTot]);
    y1_rs = reshape(y(1,f1_subss(i),:), [1 nTot]);
    ub = max([yT1_rs; y1_rs; xref(1,:)], [], 'all');
    lb = min([yT1_rs; y1_rs; xref(1,:)], [], 'all');
    range = ub - lb;
    
    subplot(2,1,i);
    plot(t, yT1_rs, t, y1_rs, t, xref(1,:), '--' );
    grid on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11); 
    ylabel(sprintf('$\\mathcal P_%i$', f1_subss(i)), 'Interpreter', 'latex');
    if i > 1
        xlabel('Time [s]', 'Interpreter', 'latex'); 
    end
    ylim([lb - 0.5*range, ub + 0.5*range]);
    legend({sprintf('$\\tilde \\theta_{%i}$', f1_subss(i)), ...
            sprintf('$\\theta_{%i}$', f1_subss(i)), ...
            sprintf('$\\theta_{%i, ref}$', f1_subss(i))} , ...
            'Interpreter', 'latex', 'Location', 'SouthWest', 'FontSize', 9);
%     linkaxes(ax,'x');
end

print('-f1', ['fig' filesep 'trajectories.eps'], '-depsc2');

%% alarms
track_bound = zeros(N, 1);
% compute alarm signals
for i=1:N
    track_bound(i) = 1.1*max(abs( squeeze(yT(1,i,t > 2 & t < tA)) - xref(1, t > 2 & t < tA)' ));
end
alarms = dcres > repmat(track_bound, [1 nTot]);
kTrans = find(t == 2);
alarms(:,1:kTrans) = zeros(N, kTrans);
dist_alarm = all(alarms(subss(nA).Ni, :), 1);
kD = find(dist_alarm, 1);
tD = t(kD);

figure(2)
set(2, 'defaultTextInterpreter', 'latex');
set(2, 'Unit', 'inches', 'Position', [0 0 6 1.9]);
% for i=1:N
%     subplot(2,3,i);
    plot(t, all(alarms(subss(nA).Ni, :), 1), 'LineWidth', 1.5, 'Color', 'red');
    grid on
    ylim([-0.1 1.1]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11); 
    yinterval = [0 1];
    yticks = linspace(yinterval(1), yinterval(2), 2); 

    ylim(yinterval + [-1, 1]*0.1);
    set(gca, 'Ytick', yticks);
    set(gca, 'Yticklabel', yticks);
    ylabel(sprintf('$\\mathcal P_%i$', nA), ...
            'Interpreter', 'latex', 'FontSize', 12);
    if i > 3
        xlabel('Time [s]', 'Interpreter', 'latex'); 
    end
% end

print('-f2', ['fig' filesep 'alarms.eps'], '-depsc2');
%% residuals 

figure(3)
set(3, 'defaultTextInterpreter', 'latex');
set(3, 'Unit', 'inches', 'Position', [0 0 12 3.5]);
for i=1:N
    rbi = max(dcres(i,t > tA), [], 2);
    subplot(2,3,i);
    plot(t, dcres(i,:).', [t(1), t(end)], repmat(track_bound(i), [1 2]), '--');
    grid on
    ylabel(sprintf('$\\|\\delta_{%i}\\| $', i), ...
            'Interpreter', 'latex', 'FontSize', 12);
    ylim([0 rbi*1.2]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11); 
    if i > 3
        xlabel('Time [s]', 'Interpreter', 'latex'); 
    end
end

print('-f3', ['fig' filesep 'residuals.eps'], '-depsc2');

%% inputs
% figure(4)
% 
% for i=1:N
%     subplot(3,2,i);
%     plot(t, reshape(u(:,i,:), [mI nTot]), ...
%          t, reshape(utilde(:,i,:), [mI nTot]), '--');
%     if i == nA
%         hold on
%         plot(t(1:end-1), [utildej{i,:}].', '--', t(1:end-1), [uj{i,:}].')
%         hold off
%     end
%     set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12); 
%     ylabel(sprintf('$\\mathcal P_%i$', i), 'Interpreter', 'latex');
% end
