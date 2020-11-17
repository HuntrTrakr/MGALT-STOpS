function [constant_value,constant_units] = constants(constant_name)
% FORM: [constant_value,constant_units] = constants(constant_name)
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -This function looks up values of various constants. The input 
% |     must be a string. The string is NOT case sensitive, but it must 
% |     be an exact match, letters-wise.
% |
% |     -To see a list of all available constants, their values, and 
% |     their units, call the function with no inputs.
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -constant_name      (1,1)       [string]        [unitless]
% |         Desired constant to look up
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -constant_value     (1,1)       [float][int]    [depends]
% |     	Numerical value of desired constant
% |     -constant_units     (1,1)       [string]        [depends]
% |       	Units of desired constant
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -References
% |         Curtis, Howard D. "Orbital Mechanics for Engineering
% |             Students", 2e --> SOI Values
% |         Wikipedia, "Standard gravitational parameter"
% |         Vallado, David A. "Fundamentals of Astrodynamics and
% |             Applications", 4e
% |
% |-----------------------------------------------------------------------



%% Celestial Bodys' Radii
names  = [ {'radius Sun'} {'radius Mercury'} {'radius Venus'} {'radius Earth'} {'radius Moon'} {'radius Mars'} {'radius Jupiter'} {'radius Saturn'} {'radius Uranus'} {'radius Neptune'} {'radius Pluto'} {'radius Eris'} ];
values = [ {695800}       {2440}             {6052}           {6378}           {1737.4}        {3390}          {69911}            {58232}           {25362}           {24622}            {1184}           {1400}          ];
units  = [ {'km'}         {'km'}             {'km'}           {'km'}           {'km'}          {'km'}          {'km'}             {'km'}            {'km'}            {'km'}             {'km'}           {'km'}          ];



%% Celestial Bodys' Sphere of Influence Diameters
names  = [ names  {'SOI Mercury'} {'SOI Venus'} {'SOI Earth'} {'SOI Moon'} {'SOI Mars'} {'SOI Jupiter'} {'SOI Saturn'} {'SOI Uranus'} {'SOI Neptune'} {'SOI Pluto'} ];
values = [ values {112000}        {616000}      {925000}      {66100}      {577000}     {48200000}      {54800000}     {51800000}     {3080000}       {112000}      ];
units  = [ units  {'km'}          {'km'}        {'km'}        {'km'}       {'km'}       {'km'}          {'km'}         {'km'}         {'km'}          {'km'}        ];



%% Celestial Bodys' Standard Gravitational Parameters
names  = [ names  {'mu Sun'}     {'mu Mercury'} {'mu Venus'} {'mu Earth'}  {'mu Moon'}  {'mu Ceres'} {'mu Mars'}  {'mu Jupiter'} {'mu Saturn'} {'mu Uranus'} {'mu Neptune'} {'mu Pluto'} {'mu Eris'}  ];
values = [ values {132712440018} {22032}        {324859}     {398600.4418} {4904.8695}  {62.6325}    {42838.37}   {126686534}    {37931187}    {5793939}     {6836529}      {871}        {1108}       ];
units  = [ units  {'km^3/s^2'}   {'km^3/s^2'}   {'km^3/s^2'} {'km^3/s^2'}  {'km^3/s^2'} {'km^3/s^2'} {'km^3/s^2'} {'km^3/s^2'}   {'km^3/s^2'}  {'km^3/s^2'}  {'km^3/s^2'}   {'km^3/s^2'} {'km^3/s^2'} ];



%% Celestial Bodys' Orbital Periods
names  = [ names  {'period Mercury'} {'period Venus'} {'period Earth'} {'period Mars'} {'period Jupiter'} {'period Saturn'} {'period Uranus'} {'period Neptune'}];
values = [ values {87.97}            {224.7}          {365.26}         {686.98}        {4332.82}          {10755.70}        {30687.15}        {60190.03}];
units  = [ units  {'days'}           {'days'}         {'days'}         {'days'}        {'days'}           {'days'}          {'days'}          {'days'}];



%% Celestial Bodys' Minimum Allowable Flyby Radius
names  = [ names  {'min rp Sun'} {'min rp Mercury'} {'min rp Venus'} {'min rp Earth'} {'min rp Mars'} {'min rp Jupiter'} {'min rp Saturn'} {'min rp Uranus'} {'min rp Neptune'} {'min rp Pluto'} ];
values = [ values {695800}       {2740}             {6351}           {6678}           {3689}          {600000}           {70000}           {40000}           {40000}            {1500}           ];
units  = [ units  {'km'}         {'km'}             {'km'}           {'km'}           {'km'}          {'km'}             {'km'}            {'km'}            {'km'}             {'km'}           ];



%% Mean Planetary Constants for Epoch J2000 (VALLADO p. 1041 Table D-3 and D-4)
names  = [ names  {'a Mercury'} {'a Venus'} {'a Earth'} {'a Mars'}  {'a Jupiter'} {'a Saturn'} {'a Uranus'} {'a Neptune'} {'a Pluto'}  ];
values = [ values {57909083}    {108208601} {149598023} {227939186} {778298361}   {1429394133} {2875038615} {4504449769}  {5915799000} ];
units  = [ units  {'km'}        {'km'}      {'km'}      {'km'}      {'km'}        {'km'}       {'km'}       {'km'}        {'km'}       ];

names  = [ names  {'e Mercury'} {'e Venus'}   {'e Earth'}   {'e Mars'}    {'e Jupiter'} {'e Saturn'}  {'e Uranus'}  {'e Neptune'} {'e Pluto'} ];
values = [ values {0.205631752} {0.006771882} {0.016708617} {0.093400620} {0.048494851} {0.055508622} {0.046295898} {0.008988095} {0.249050} ];
units  = [ units  {'none'}      {'none'}      {'none'}      {'none'}      {'none'}      {'none'}      {'none'}      {'none'}      {'none'} ];

names  = [ names  {'i Mercury'} {'i Venus'} {'i Earth'} {'i Mars'}   {'i Jupiter'} {'i Saturn'} {'i Uranus'} {'i Neptune'} {'i Pluto'}   ];
values = [ values {7.00498625}  {3.9446619} {0}         {1.84972648} {1.30326966}  {2.48887810} {0.77319617} {1.76995221}  {17.14216667} ];
units  = [ units  {'°'}         {'°'}       {'°'}       {'°'}        {'°'}         {'°'}        {'°'}        {'°'}         {'°'}         ];

names  = [ names  {'RAAN Mercury'} {'RAAN Venus'} {'RAAN Earth'} {'RAAN Mars'} {'RAAN Jupiter'} {'RAAN Saturn'} {'RAAN Uranus'} {'RAAN Neptune'} {'RAAN Pluto'} ];
values = [ values {48.33089304}    {76.67992019}  {0}            {49.55809321} {100.46444064}   {113.66552370}  {74.00594723}   {131.78405702}   {110.29713889} ];
units  = [ units  {'°'}            {'°'}          {'°'}          {'°'}         {'°'}            {'°'}           {'°'}           {'°'}            {'°'}          ];

names  = [ names  {'omega Mercury'} {'omega Venus'} {'omega Earth'} {'omega Mars'} {'omega Jupiter'} {'omega Saturn'} {'omega Uranus'} {'omega Neptune'} {'omega Pluto'} ];
values = [ values {77.45611904}     {131.56370724}  {102.93734808}  {336.06023398} {14.33130924}     {93.05678728}    {173.00515922}   {48.12369050}     {224.13486111}  ];
units  = [ units  {'°'}             {'°'}           {'°'}           {'°'}          {'°'}             {'°'}            {'°'}            {'°'}             {'°'}           ];

names  = [ names  {'lambda_M Mercury'} {'lambda_M Venus'} {'lambda_M Earth'} {'lambda_M Mars'} {'lambda_M Jupiter'} {'lambda_M Saturn'} {'lambda_M Uranus'} {'lambda_M Neptune'} {'lambda_M Pluto'} ];
values = [ values {252.25090551}       {181.97980084}     {100.46644851}     {355.43327463}    {34.35148392}        {50.07747138}       {314.05500511}      {304.34866548}       {238.74394444}     ];
units  = [ units  {'°'}                {'°'}              {'°'}              {'°'}             {'°'}                {'°'}               {'°'}               {'°'}                {'°'}              ];



%% Miscellaneous Earth Information
names  = [ names  {'eccentricity of Earth'} ];
values = [ values {0.081819221456}          ];
units  = [ units  {'none'}                  ];



%% Astronomical Units
names  = [ names  {'AU'}     {'TU'}   ];
values = [ values {1.496e8}  {58.13}  ];
units  = [ units  {'km'}     {'days'} ];



%% If No Inputs Show Available Inputs
if nargin < 1
    clc
    disp([names',values',units']);
    return
end



%% Else, Retrieve Value
index = find(strcmpi(constant_name,names));
% index = strcmpi(constant_name,names);       % Should slightly speed up the code

if index
    constant_value = values{index};
    constant_units = units {index};
else
    errorPathDisplay();
    fprintf(2,'Constant not found. Try a different name.\n')
    return
end



end


