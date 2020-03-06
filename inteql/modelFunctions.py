import numpy as np
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.dummy import DummyRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from scipy.stats import linregress


def decision_tree_regressor(data, X_label, y_label, random_state, test_size=0.3, max_depth=None):
    """
    This function generate a decision tree regressor. 
    First it splits the data into train and test. 
    Then, treins the model.
    Fnally, predicts the y values for the test set and compute its error.
    It receives:
        - Data: as a data frame
        - List of the feature column names: list
        - Target column name: str
        - Random state: float
        - The size of the test set: float [0,1]
        - The maximum depth of the tree: int
    It returns a dictionary:
        - 'rmse': The root mean square error of the predictor
        - 'model': The trained model
        - 'y_pred': The predictions of the test set
    """
    
    # Split data
    X = data[X_label]
    y = data[y_label]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    
    # print(X_train.shape)
    # Create a model and train it
    model = DecisionTreeRegressor(random_state=random_state, max_depth=max_depth)
    model = model.fit(X_train, y_train)
    
    # Predict values for test set and asses error
    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    
    # Return result object
    return {'rmse':rmse, 'model':model, 'y_pred':y_pred}

def random_forest_regressor(data, X_label, y_label, random_state, test_size=0.3, max_depth=None, n_estimators=100):
    """
    This function generate a Random Forest Regressor. 
    First it splits the data into train and test. 
    Then, treins the model.
    Fnally, predicts the y values for the test set and compute its error.
    It also computes the regression line for the test and predictions.
    It receives:
        - Data: as a data frame
        - List of the feature column names: list
        - Target column name: str
        - Random state: float
        - The size of the test set: float [0,1]
        - The maximum depth of the tree: int
        - The number of estimators: int. By default is 100.
    It returns a dictionary:
        - 'rmse': The root mean square error of the predictor
        - 'model': The trained model
        - 'y_pred': The predictions of the test set
        - 'r_value': The R value of the regression line
        - 'y_test': The target column of the test set
    """
    
    # Split data
    X = data[X_label]
    y = data[y_label]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    
    # Create a model and train it
    model = RandomForestRegressor(random_state=random_state, max_depth=max_depth, n_estimators=n_estimators)
    model = model.fit(X_train, y_train)
    
        
    # Predict values for test set and asses error
    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    slope, intercept, r_value, p_value, std_err = linregress(y_test, y_pred)
    
    
    # Return result object
    return {'rmse':rmse, 'model':model, 'y_pred':y_pred, 'r_value':r_value, 'y_test':y_test}

def dummy_regressor(data, X_label, y_label, random_state, test_size=0.3):
    """
    This function generate a Dummy regressor. 
    First it splits the data into train and test. 
    Then, trains the model.
    Finally, predicts the y values for the test set and compute its error.
    Then, trains the model.
    Finally, predicts the y values for the test set and compute its error.
    It receives:
        - Data: as a data frame
        - List of the feature column names: list
        - Target column name: str
        - Random state: float
        - The size of the test set: float [0,1]
        - The maximum depth of the tree: int
    It returns a dictionary:
        - 'rmse': The root mean square error of the predictor
        - 'model': The trained model
        - 'y_pred': The predictions of the test set
    """
    
    # Split data
    X = data[X_label]
    y = data[y_label]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    
    # Create a model and train it
    model = DummyRegressor()
    model = model.fit(X_train, y_train)
    
    # Predict values for test set and asses error
    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))

    # Return result object
    return {'rmse': rmse, 'model': model, 'y_pred': y_pred}':y_test}