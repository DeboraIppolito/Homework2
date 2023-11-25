#include "kdl_ros_control/kdl_control.h"

KDLController::KDLController(KDLRobot &_robot)
{
    robot_ = &_robot;
}

//CONTROLLO JOINT SPACE
Eigen::VectorXd KDLController::idCntr(KDL::JntArray &_qd,
                                      KDL::JntArray &_dqd,
                                      KDL::JntArray &_ddqd,
                                      double _Kp, double _Kd)
{
    // read current joint state
    Eigen::VectorXd q = robot_->getJntValues();
    Eigen::VectorXd dq = robot_->getJntVelocities();

    // calculate errors
    Eigen::VectorXd e = _qd.data - q;
    Eigen::VectorXd de = _dqd.data - dq;

    Eigen::VectorXd ddqd = _ddqd.data;
    return robot_->getJsim() * (ddqd + _Kd*de + _Kp*e)
            + robot_->getCoriolis() + robot_->getGravity() /*friction compensation?*/;
}


//CONTROLLO OPERATIONAL SPACE
Eigen::VectorXd KDLController::idCntr(KDL::Frame &_desPos,
                                      KDL::Twist &_desVel,
                                      KDL::Twist &_desAcc,
                                      double _Kpp, double _Kpo,
                                      double _Kdp, double _Kdo)
{
   
   // calculate gain matrices
   Eigen::Matrix<double,6,6> Kp, Kd;
   Kp.block(0,0,3,3) = _Kpp*Eigen::Matrix3d::Identity(); //costruisco la matrice 3x3 dei guadagni sull'errore di posizione
   Kp.block(3,3,3,3) = _Kpo*Eigen::Matrix3d::Identity(); //costruisco la matrice 3x3 dei guadagni sull'errore di orientamento
   Kd.block(0,0,3,3) = _Kdp*Eigen::Matrix3d::Identity();//costruisco la matrice 3x3 dei guadagni sulla derivata dell'errore di posizione
   Kd.block(3,3,3,3) = _Kdo*Eigen::Matrix3d::Identity();//costruisco la matrice 3x3 dei guadagni sulla derivata dell'errore di orientamento


   //conversione da KDL::Jacobian a Eigen::Matrix e calcolo della pseudoinversa
   Eigen::Matrix<double,6,7> J = robot_->getEEJacobian().data;
   Eigen::Matrix<double,7,6> Jpinv = pseudoinverse(J);

   Eigen::Matrix<double,7,7> I = Eigen::Matrix<double,7,7>::Identity();
   Eigen::Matrix<double,7,7> M = robot_->getJsim();
   //Eigen::Matrix<double,7,6> Jpinv = weightedPseudoInverse(M,J);


   // posizione desiderata ed effettiva
   Eigen::Vector3d p_d(_desPos.p.data);
   Eigen::Vector3d p_e(robot_->getEEFrame().p.data);
   Eigen::Matrix<double,3,3,Eigen::RowMajor> R_d(_desPos.M.data);
   Eigen::Matrix<double,3,3,Eigen::RowMajor> R_e(robot_->getEEFrame().M.data);
   R_d = matrixOrthonormalization(R_d);
   R_e = matrixOrthonormalization(R_e);

   // velocita' desiderata ed effettiva
   Eigen::Vector3d dot_p_d(_desVel.vel.data);
   Eigen::Vector3d dot_p_e(robot_->getEEVelocity().vel.data);
   Eigen::Vector3d omega_d(_desVel.rot.data);
   Eigen::Vector3d omega_e(robot_->getEEVelocity().rot.data);

   // accelerazione desiderata ed effettiva
   Eigen::Matrix<double,6,1> dot_dot_x_d;
   Eigen::Matrix<double,3,1> dot_dot_p_d(_desAcc.vel.data);
   Eigen::Matrix<double,3,1> dot_dot_r_d(_desAcc.rot.data);

   //calcolo errore di posizione e velocità
   Eigen::Matrix<double,3,1> e_p = computeLinearError(p_d,p_e);
   Eigen::Matrix<double,3,1> dot_e_p = computeLinearError(dot_p_d,dot_p_e);   
   //calcolo errore orientamento 
   Eigen::Matrix<double,3,1> e_o = computeOrientationError(R_d,R_e);
   Eigen::Matrix<double,3,1> dot_e_o = computeOrientationVelocityError(omega_d, omega_e, R_d, R_e);

   //ERRORE COMPLESSIVO                                                        
   Eigen::Matrix<double,6,1> x_tilde;
   Eigen::Matrix<double,6,1> dot_x_tilde;
   x_tilde << e_p, e_o;
   dot_x_tilde << dot_e_p, -omega_e;
   dot_dot_x_d << dot_dot_p_d, dot_dot_r_d;

   // null space control
   double cost;
   Eigen::VectorXd grad = gradientJointLimits(robot_->getJntValues(),robot_->getJntLimits(),cost);


    //matrice contenente y = xd_dot_dot - J_dot*q_dot + Kd*x_tilde_dot + Kp*x_tilde
    Eigen::Matrix<double,6,1> y;
    y << dot_dot_x_d - robot_->getEEJacDotqDot() + Kd*dot_x_tilde + Kp*x_tilde;

    //restituiamo l'ingresso di controllo u = By + n
       return M * (Jpinv*y+ (I-Jpinv*J)*(/*- 10*grad */- 1*robot_->getJntVelocities()))+ robot_->getGravity() + robot_->getCoriolis();
}

