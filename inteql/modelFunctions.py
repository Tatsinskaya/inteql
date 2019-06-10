import numpy as np
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.dummy import DummyRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from scipy.stats import linregress


def decision_tree_regressor(data, X_label, y_label, random_state, test_size=0.3, max_depth=None):
    """
    This bla bla bla
     clsd
    """
    
    # Split data
    X = data[X_label]
    y = data[y_label]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    
    print(X_train.shape)
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
    This bla bla bla
    
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

def dummy_regressor(data, y_label, random_state, test_size=0.3):
    """
    This bla bla bla
    
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
    return {'rmse':rmse, 'model':model, 'y_pred':y_pred}