function td = calc_td(tp,DC)
% calculates the interpuls delay 'td'
% input: puls duration 'tp', duty cycle 'DC'
% last change: 2014/04/02 by PS

td=tp*(1/DC-1);
