import pandas as pd
import seaborn as sns



# see notebook for full set of code - this is a test

# Load the data
data = pd.read_csv('https://raw.githubusercontent.com/FBosler/AdvancedPlotting/master/combined_set.csv')
# this assigns labels per year
data['Mean Log GDP per capita']  = data.groupby('Year')['Log GDP per capita'].transform(
    pd.qcut,
    q=5,
    labels=(['Lowest','Low','Medium','High','Highest'])
)

#Plotting with seaborn
sns.reset_defaults()
sns.set(
    rc={'figure.figsize':(7,5)},
    style="white" # nicer layout
)

#Single-variable distribution
sns_data=data[
    (data['Year']==2018) &
    (data['Continent']=='Asia')
]

sns_plot=sns.distplot(
    sns_data['Life Ladder'],
    label='Life Ladder'
)


sns_plot.figure.savefig("output.png")