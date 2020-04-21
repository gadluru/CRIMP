function shifts_interp = interp_rigid_shifts(shifts,respiratory_signal,reference_points)

shifts_interp = zeros(size(shifts,1),length(respiratory_signal));
shift_ref = respiratory_signal(reference_points);

shifts = squeeze(sum(shifts,2));
for i=1:size(shifts,1)
    f = fit(double(shift_ref), double(shifts(i,:)'), 'poly1' );
    shifts_interp(i,:) = f(respiratory_signal);
end