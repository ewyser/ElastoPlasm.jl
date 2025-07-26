## Deformation \& Strain
The theoretical background and general notions related either to Continuum Mechanics, Infinitesimal and Finite Deformation Mechanics or Elasticity and Elasto-plasticity theories presented in this chapter can be found in the following classical text book references: \cite{Bower2009,Hashiguchi2012,Spencer2012,Borja2013,Doghri2013}. The purpose is to avoid too many references within the text body and to improve the reading.
\section{Motions and deformations}
\subsection{Description of motion}
Let us consider a material point in a continuum body defined by its coordinates $\boldsymbol{X}$ in a reference configuration $\Omega_0 \subset \mathbb{R}^3$ at a reference time $t_0=0$. For a given later time $t>0$ the material point coordinate $\boldsymbol{X}$ has now coordinates $\boldsymbol{x}$ in the current configuration $\Omega_t \subset \mathbb{R}^3$, then the equation
\begin{align}
	\boldsymbol{x}=\boldsymbol{\phi}(\boldsymbol{X},t), \quad \text{\eg} \quad x_i = \phi_i(X_I,t),
\end{align}
describes a \textit{motion} of the continuum where $\boldsymbol{\phi}$ is a transformation (\eg a mapping) from the reference configuration to the current configuration. 
\begin{figure}[htbp]
	\centering
	\includegraphics[scale=0.6]{figs/motion_description}
	\caption{Coordinate $\boldsymbol{X}$ of a material point in the reference configuration $\Omega_0$ and its updated coordinate $\boldsymbol{x}=\boldsymbol{\phi}(\boldsymbol{X},t)$ in the current configuration, in a fixed Cartesian orthonormal frame $(\mathcal{O},\boldsymbol{e}_1,\boldsymbol{e}_2,\boldsymbol{e}_3)$.}
	\label{motion}
\end{figure}

Any given material point is specified by its position vectors, \eg 
\begin{align}
	\boldsymbol{X} &= \{X_I\}, \quad &\text{in the \textit{reference} configuration} \\ 
	\boldsymbol{x} &= \{x_i\}, \quad  &\text{in the \textit{current} configuration} 
\end{align} 
and upper cases (\eg $I,J,K$) and lower cases (\eg $i,j,k$) indices designate the coordinates in the reference and current configurations, respectively. By convention, $X_I$ are called \textit{the material coordinates} whereas $x_i$ are called \textit{the spatial coordinates}\footnote{Similarly, the material coordinates can be referred to as \textit{Lagrangian coordinates} and, the spatial coordinates can be referred to as \textit{Eulerian coordinates}}.

During a rigid-body motion made of both translation and rotation, a body moves without changing its shape; the distance and the orientation between two material points do not change during the motion. A \textit{translation} is a rigid-body motion during which every material point undergoes the same displacement. Such motion is described by the following equation
\begin{align}
	\boldsymbol{x}=\boldsymbol{X}+\boldsymbol{c}(t),
\end{align}
where the vector $\boldsymbol{c}(t)$ is frame invariant and only depends on $t$. A \textit{rotation} is also a rigid-body motion during which every material point undergoes the same rotation about any arbitrary axis of origin $\mathcal{O}$. It is described by the following equation
\begin{align}
	\boldsymbol{x}=\boldsymbol{Q}(t) \boldsymbol{X},
	\label{rigid_body_motion_eq}
\end{align}
where $\boldsymbol{Q}(t)$ is an orthogonal tensor. 

As such, any rigid-body motion is a combination of both a translational part and a rotation part about an axis and, it can be described by equations of the form 
\begin{align}
	\boldsymbol{x}&=\boldsymbol{Q}(t) \boldsymbol{X} + \boldsymbol{c}(t), \quad \text{or} \quad \boldsymbol{X}=\boldsymbol{Q}^T(t) \boldsymbol{x} - \boldsymbol{Q}^T(t) \boldsymbol{c}(t).
\end{align} 

\subsection{Deformation}
A body will change its shape as well as its position (translation) and its orientation (rotation) during any motion. A motion during which a change in shape takes place is called a deformation, independently of a change of position or orientation. 

As such, the main problem in deformation analysis is to separate rigid-body motion from deformation, which has to be \textit{invariant} with respect to rigid-body motion.

The \textit{deformation gradient} is defined by:
\begin{align}
	\boldsymbol{F} \equiv \dfrac{\partial \boldsymbol{\phi}}{\partial \boldsymbol{X}} = \dfrac{\partial \boldsymbol{x}}{\partial \boldsymbol{X}}, \quad \text{\eg} \quad F_{iJ}=\dfrac{\partial x_i}{\partial X_J},
\end{align}
and describes how much the current coordinate $x_i$ varies w.r.t to the reference coordinate $X_J$. Taking the inverse of the deformation gradient, \eg $\boldsymbol{F}^{-1}$ gives the inverse relation between the reference coordinate w.r.t the current configuration. If there is no motion, then $x_i=X_I$ and so $F_{iI}$ reduces to $\delta_{iI}$, \eg $F_{iI}=\delta_{iI}$. 

Alternatively, the current position is defined by $\boldsymbol{x}=\boldsymbol{X}+\boldsymbol{u}$, where $\boldsymbol{u}$ is the displacement. By definition, the deformation gradient is also given by
\begin{align}
	\begin{aligned}
		\boldsymbol{F} &= \dfrac{\partial}{\partial \boldsymbol{X}}(\boldsymbol{X}+\boldsymbol{u}),\\
		\boldsymbol{F} &= \dfrac{\partial \boldsymbol{X}}{\partial \boldsymbol{X}}+\dfrac{\partial \boldsymbol{u}}{\partial \boldsymbol{X}},\\
		\boldsymbol{F} &= \boldsymbol{I}+\dfrac{\partial \boldsymbol{u}}{\partial \boldsymbol{X}} \label{dF_def},
	\end{aligned}
\end{align}
and rearranging terms, the \textit{displacement} gradient tensor is given by
\begin{align}
	\dfrac{\partial \boldsymbol{u}}{\partial \boldsymbol{X}} = \boldsymbol{F}-\boldsymbol{I},  \quad \text{\eg} \quad \dfrac{\partial u_i}{\partial X_J}=\dfrac{\partial x_i}{\partial X_J}-\delta_{iJ}.
\end{align}

However, the displacement gradient is expressed w.r.t the reference configuration. Differentiating the displacement w.r.t to the current configuration yields
\begin{align}
	\dfrac{\partial \boldsymbol{u}}{\partial \boldsymbol{x}} = \boldsymbol{I} - \boldsymbol{F}^{-1},  \quad \text{\eg} \quad \dfrac{\partial u_i}{\partial x_j}=\delta_{iJ}-\dfrac{\partial X_I}{\partial x_j}.
\end{align}

Let consider an infinitesimal line element vector $d\boldsymbol{X}$ with origin $\boldsymbol{X}$ in the reference configuration and an infinitesimal line element vector $d\boldsymbol{x}$ with origin $\boldsymbol{x}$ in the current configuration (see Fig. \ref{deformation_description}). By definition, one has
\begin{align}
	d\boldsymbol{x} = \boldsymbol{F} d\boldsymbol{X} \quad \text{\eg} \quad dx_i = F_{iJ} dX_J,
\end{align}  
and demonstrates that the deformation gradient is the fundamental measure of deformation in continuum mechanics. 

\begin{figure}[htbp]
	\centering
	\includegraphics[scale=0.6]{figs/deformation_description}
	\caption{Deformation of an infinitesimal line element $d\boldsymbol{X}$ in a reference configuration.}
	\label{deformation_description}
\end{figure}

\subsection{Finite deformation and strain measure}
A \textit{measure} of the deformation (\eg the \textit{strain}) should be \textit{invariant by rotation}. The deformation gradient tensor $\boldsymbol{F}$ does not have to satisfy this property. In fact, in the rigid-body motion given by Eq. \ref{rigid_body_motion_eq}, the deformation gradient reduces to $\boldsymbol{F}=\boldsymbol{Q}(t)$. As such, $\boldsymbol{F}$ itself is not a suitable measure of the deformation under a rigid-body motion and, it does not allow a proper invariant definition of the strain measure. 

What is then a suitable measure of deformation ? The answer lies in the definition of additional deformation tensors which are briefly presented in the following.

Using the polar decomposition theorem, the deformation gradient can be decomposed, in a multiplicative manner, by a product of two second-order tensors. It follows that 
\begin{align}
	\boldsymbol{F}=\boldsymbol{R} \boldsymbol{U} = \boldsymbol{V} \boldsymbol{R} \quad \text{\eg} \quad F_{iJ} = R_{iK}U_{KJ} = V_{ik}R_{kJ},
\end{align}  
where $\boldsymbol{R}$ is a proper orthogonal tensor (\eg a rotation tensor), and $\boldsymbol{U}$ and $\boldsymbol{V}$ are the right stretch (\eg defined w.r.t. the reference configuration in the \textit{material} coordinates) and left stretch (\eg defined w.r.t the current configuration in the \textit{spatial} coordinates) tensors, respectively. The polar decomposition \textit{multiplicatively} decomposes the deformation gradient into orthogonal (\eg rotation) and stretch tensors.

The right Cauchy-Green deformation tensor $\boldsymbol{C}$ is a first suitable measure of deformation and is given by 
\begin{align}
	\boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F}= \boldsymbol{U}^T \boldsymbol{R}^T \boldsymbol{R} \boldsymbol{U} = \boldsymbol{U}^2 \quad \text{\eg} \quad C_{IJ}=\dfrac{\partial x_k}{\partial X_I} \dfrac{\partial x_k}{\partial X_J}.
\end{align}
and, since $\boldsymbol{R}$ is orthogonal, hence $\boldsymbol{R}\boldsymbol{R}^T=\boldsymbol{R}^T\boldsymbol{R}=\boldsymbol{I}$, which demonstrates the invariance by rotation of $\boldsymbol{C}$.

Similarly, the left Cauchy-Green deformation tensor $\boldsymbol{b}$ is a second suitable measure of deformation and is given by 
\begin{align}
	\boldsymbol{b} = \boldsymbol{F}  \boldsymbol{F}^T = \boldsymbol{V} \boldsymbol{R} \boldsymbol{R}^T \boldsymbol{V}^T = \boldsymbol{V}^2 \quad \text{\eg} \quad b_{ij}=\dfrac{\partial x_i}{\partial X_K} \dfrac{\partial x_j}{\partial X_K}.
\end{align}

These two tensors allow to define the \textit{Lagrangian strain} tensor $\boldsymbol{\gamma}$ (\eg the strain is defined w.r.t the reference configuration) and the \textit{Eulerian strain} tensor $\boldsymbol{\eta}$ (\eg the strain is defined w.r.t the current configuration), respectively, and are defined as
\begin{align}
	\boldsymbol{\gamma}&=\dfrac{1}{2}(\boldsymbol{C}-\boldsymbol{I}),\\
	\boldsymbol{\eta}&=\dfrac{1}{2}(\boldsymbol{I}-\boldsymbol{b}^{-1}).
\end{align}

The expressions of both $\boldsymbol{\gamma}$ and  $\boldsymbol{\eta}$ can be expressed in terms of the displacement gradient, which results in the following
\begin{align}
	\gamma_{IJ}&=\dfrac{1}{2}\left( \dfrac{\partial u_I}{\partial X_J}+\dfrac{\partial u_J}{\partial X_I}+\dfrac{\partial u_k}{\partial X_I}\dfrac{\partial u_k}{\partial X_J} \right), \label{lagrangian_tensor}\\
	\eta_{ij}&=\dfrac{1}{2}\left( \dfrac{\partial u_i}{\partial x_j}+\dfrac{\partial u_j}{\partial x_i}-\dfrac{\partial u_K}{\partial x_i}\dfrac{\partial u_K}{\partial x_j} \right).
	\label{eulerian_tensor}
\end{align}

\subsection{Infinitesimal strain}
The fundamental assumption of the infinitesimal strain theory is that the current configuration is not significantly different from the reference configuration. This implies that all the components of the displacement gradient are numerically small, \eg $| \partial u_i / \partial X_J |  \ll 1$ and, ii) the squares and products of these quantities are neglected. 

To demonstrate the equivalence \citep{Spencer2012} , let consider the displacement gradient expressed in the spatial coordinates, \eg
\begin{align}
	\dfrac{\partial \boldsymbol{u}}{\partial\boldsymbol{x}}=\boldsymbol{I}-\boldsymbol{F}^{-1} \quad \text{\eg} \quad \dfrac{\partial u_i}{\partial x_j} = \delta_{ij} - \dfrac{\partial X_I}{\partial x_j}, 
\end{align}
where, by binomial expansion, $\boldsymbol{I}-\boldsymbol{F}^{-1}$ yields to
\begin{align}
	\boldsymbol{I}-\boldsymbol{F}^{-1} = \boldsymbol{I} - \left( \boldsymbol{I}-(\boldsymbol{F}-\boldsymbol{I}) + (\boldsymbol{F}-\boldsymbol{I})^2 -(\boldsymbol{F}-\boldsymbol{I})^3+...\right).
\end{align}
The displacement gradient now reads 
\begin{align}
	\dfrac{\partial \boldsymbol{u}}{\partial\boldsymbol{x}}= (\boldsymbol{F}-\boldsymbol{I}) - (\boldsymbol{F}-\boldsymbol{I})^2 + (\boldsymbol{F}-\boldsymbol{I})^3-..., 
\end{align}
and, considering $\boldsymbol{F}-\boldsymbol{I} = \partial \boldsymbol{u}/\partial \boldsymbol{X}$ and neglecting higher order terms, the displacement gradient formulated in the material coordinates reduces to
\begin{align}
	\dfrac{\partial \boldsymbol{u}}{\partial\boldsymbol{x}} = \dfrac{\partial \boldsymbol{u}}{\partial\boldsymbol{X}} \quad \text{\eg} \quad \dfrac{\partial u_i}{\partial x_j} = \dfrac{\partial u_i}{\partial X_J}.
\end{align}

To first order in the displacement gradient, it follows from Eqs. \ref{lagrangian_tensor}-\ref{eulerian_tensor} that $\gamma_{IJ} \simeq \eta_{ij}$. The \textit{infinitesimal strain} tensor $\boldsymbol{\epsilon}$ is defined as 
\begin{align}
	\boldsymbol{\epsilon} = \dfrac{1}{2} \left( \boldsymbol{F} + \boldsymbol{F}^T \right)- \boldsymbol{I} \quad \text{\eg} \quad \epsilon_{ij}=\dfrac{1}{2} \left(\dfrac{\partial u_i}{\partial X_J}+ \dfrac{\partial u_j}{\partial X_I} \right). \label{infinitesiaml_strain_df}
\end{align}
Hence, the infinitesimal strain tensor can be regarded as an exact but first order formulation of the displacement gradient tensor. 
