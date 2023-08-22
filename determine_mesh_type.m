function canal_flag = determine_mesh_type

[selectedButton,dlgShown]=uigetpref('mygraphics',... % Group
       'get_mesh_type',...           % Preference
       'User Selection',...                    % Window title
       {'Choose Mesh Type '
        ''
        'Typically, solid meshes are used for epiphyseal sections,'
        'and donut meshes are used for midshaft sections'},...
       {'solid','donut';'solid','donut'},...       % Values and button strings
       'DefaultButton','donut','CheckboxString','') ;          % Default choice
if selectedButton=='donut'
    canal_flag=1;
else
    canal_flag=0;
end
